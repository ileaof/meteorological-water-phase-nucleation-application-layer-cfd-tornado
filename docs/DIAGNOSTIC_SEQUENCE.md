# Instrumentação passiva e diagnóstico temporal 3D

O objetivo é diagnosticar a cadeia dinâmica, sem alterar física, solver numérico ou termodinâmica. O núcleo recebe somente chamadas opcionais a um observador. Sem `sim.diagnostic_observer`, elas não capturam dados. O observador fornecido lê cópias dos campos; os testes verificam igualdade bit a bit dos campos com e sem observador, em CPU/GPU e no caminho de baixa memória.

## Executar e reproduzir

Dependência opcional: `h5py` (`pip install -e .[diagnostics]`). NumPy, SciPy e Matplotlib já pertencem ao projeto. A execução GPU utiliza a instalação CuPy existente; não ocorre fallback silencioso para CPU.

```powershell
python scripts/run_diagnostic_sequence.py --out outputs/novo_diagnostico --duration 3900 --interval 30 --device gpu
python scripts/analyze_diagnostic_sequence.py outputs/novo_diagnostico/sequence.h5
```

O diretório de execução precisa ser novo. A captura HDF5 usa criação exclusiva e não sobrescreve resultados. Os instantes de saída são o primeiro fim de passo após o intervalo solicitado: a instrumentação não altera dt para acertar a cadência. Todos os horários reais e passos estão registrados. O último instante inclui o passo final já previsto pela duração da integração. Se faltar espaço ou surgir uma falha, o status é `interrupted`, que o analisador rejeita. A análise só abre o arquivo depois de fechado.

O caso padrão é o pai idealizado de `scratchpad/tornado_intensity_L_gpu.py`: domínio 72 × 72 × 15 km, 120 × 120 × 48, estiramento vertical 1.05, hodógrafa quarter-circle com U_max=30 m/s, bolha de 5 K, LES C_s=0.20, arrasto existente com divergência de tensão nos primeiros 150 m e lei logarítmica com z0=0.1 m. Há transformação para o referencial da tempestade por Bunkers, como no produtor original. **Não há ajuste desses parâmetros em função dos resultados.** A condição inicial é gerada do zero. Nenhum cache sem metadados é usado.

## Esquema do arquivo `sequence.h5`

| Grupo | Conteúdo e semântica |
|---|---|
| Atributos raiz | Versão do esquema, status, configuração, fontes e hashes, comando, plataforma, equipamento, semântica das duas pressões |
| `grid` | xc/yc/zc e xf/yf/zf em metros; origem fixa; periodicidade |
| `base` | Perfis do estado base efetivamente usado: theta0, qv0, p0, T0, rho0, u0, v0, zc; rho0 nas faces verticais; f |
| `steps` | t inicial, dt e índice de **cada** passo, sem subamostragem |
| `snapshots/NNNNN` | Estado no fim do passo, tempo real, início do intervalo e atributo complete |
| Campos do snapshot | u/v/w **nativos nas faces**, p, p_dyn se presente, theta, qv, ql, qi, qr, qs, qg, qh, T, rho, RH_w, na coluna inteira |
| `increments/etapa/u,v,w` | Soma das mudanças nativas de velocidade aplicadas pela etapa desde a saída anterior, em m/s; dividir pelo intervalo real dá tendência média em m/s² |
| `kinematic_integrals` | Integrais em dt dos termos vetoriais de advecção, alongamento e dilatação e dos termos verticais separados de estiramento e inclinação, em s^-1 |

Os incrementos e integrais cobrem todos os x/y e os primeiros 2 km, mais três níveis de halo vertical. Os estados, velocidades e pressões cobrem os 15 km, para não impedir o rastreamento de parcelas pela coluna. Todas as quantidades salvas preservam a precisão dos arrays; compressão gzip nível 1 e shuffle são sem perda. Os termos integrados são **médias de intervalo**, não tendências instantâneas de seu último passo. As pressões representam o último estado projetado, e as variáveis termodinâmicas o fim do transporte/microfísica.

Etapas registradas: contornos iniciais; LES (incluindo retirada/restauração do vento base); advecção; flutuabilidade; Coriolis; arrasto superficial; forçantes externas; guardas; contornos do preditor; projeção; contornos pós-projeção; transporte/microfísica/contornos. Esta última etapa registra qualquer mudança de velocidade feita nesse trecho, e não alega que o transporte escalar gere diretamente uma tendência de momento. Mudanças termodinâmicas ficam nos snapshots e afetam a flutuabilidade nos passos posteriores.

Esta primeira implementação destina-se ao **pai de malha fixa**. Não se deve anexá-la diretamente a uma hierarquia móvel/composta e afirmar sincronização multínivel. Alterações de velocidade fora das etapas observadas são detectadas na abertura do passo. Origem, contornos e histórico de hierarquia precisariam de um observador específico para serem diagnosticados como tal; essa execução não contém ninhos.

## Orçamento independente

O analisador não chama os operadores de derivação ou classificação do solver. Usa médias das faces para os centros e diferenças de segunda ordem de NumPy nas coordenadas reais, incluindo z não uniforme. As estatísticas excluem duas células laterais e os halos acima de 2 km. O residual nativo de massa usa diferenças nas faces, sem recentramento.

Para omega=(xi,eta,zeta)=curl(u), a forma material relativa é:

```text
D omega/Dt = (omega · grad)u − omega div(u) + curl(forças)
S_z = zeta dw/dz
T_z = xi dw/dx + eta dw/dy
D_z = −zeta div(u)
A = −(u · grad)omega
```

As contribuições planetárias não são escondidas em omega: entram pelo rotacional do incremento Coriolis realmente aplicado. Os termos contínuos são amostrados antes da advecção em todos os passos e ponderados por dt. O curl dos incrementos nativos por operador é calculado com o mesmo operador linear independente do início ao fim do intervalo.

Há dois fechamentos deliberadamente separados:

1. **Discreto:** curl(u_final) − curl(u_inicial) menos a soma dos curls dos incrementos. Deve fechar a arredondamento; isso verifica a contabilidade, não prova validade física.
2. **Contínuo/material:** termos S, T e D, rotacionais das outras etapas e o restante da advecção comparados à mudança ao longo das parcelas. O restante da advecção é curl(Δu_advecção) − integral(A+alongamento+dilatação). Inclui diferenças entre formas conservativa e material em escoamento divergente, discretização, limitadores e divisão de etapas. Não é rotulado como difusão numérica pura.

As forças de flutuabilidade produzem vorticidade horizontal por seu rotacional. Seus efeitos sobre zeta via inclinação são examinados com a trajetória completa de xi/eta/zeta, em vez de usar um termo vertical baroclínico nulo para declarar ausência de baroclinicidade. Não foi feita decomposição inversa da pressão nem alteração de sua participação na termodinâmica.

## Identidade, trajetórias e incerteza

A semente é predefinida: máximo de zeta × max(w,0) entre 500 e 1500 m perto de t=1800 s, pertencente a um componente de zeta ≥0.003 s^-1 conectado por seis vizinhos. O componente é associado para frente e para trás por sobreposição de voxels. Exige pelo menos três células compartilhadas e IoU≥0.02; competidores com ≥80% do melhor escore geram `ambiguous` e encerram o rastreamento naquela direção. Perda não dispara nova seleção do máximo do domínio. A regra fornece continuidade Euleriana, não a identidade material exata em fusões. O limiar é documentado e não classifica tornado.

As sementes lagrangianas são uma grade determinística em torno do pico baixo do componente. O integrador RK4 usa as faces nativas, interpolação espacial trilinear e temporal linear. Executa comparações de subpasso 2 s/1 s e de cadência de saída integral/alternada. Não extrapola abaixo do primeiro centro vertical ou fora do domínio. Mesmo em domínio periódico, as trajetórias são censuradas na borda para evitar interpretar um retorno artificial como origem comprovada.

O orçamento ao longo das parcelas interpola as integrais de intervalo no ponto médio da trajetória entre saídas. Seu residual representa também erro de amostragem e quadratura; não deve ser confundido com o fechamento discreto. A análise conserva os três componentes, assinala intervalos fora de 0–2 km e limita as origens/contribuições de referência ao trecho válido **contíguo** que contém a semente. Os dados completos permitem auditar outros trechos sem unir lacunas.

`p` termodinâmico e `p_dyn` da projeção ficam separados. O condensado visível usa somente ql+qi, sem precipitação somada. Um disco de 1.2 km e anel de 3–6 km definem a anomalia de pressão na altura do pico; isso não isola exclusivamente pressão de um tornado. A ausência de condensado em nível baixo e a existência de depressão dinâmica são reportadas independentemente.

## Entrega e testes

O analisador gera `analysis/RELATORIO_FINAL.md`, com categorias **comprovado**, **suportado** e **sugestivo/não identificado**, além de JSON, CSV e figuras de orçamento, componente rastreado, coerência vertical, trajetórias e pressão versus condensado. Não há recomendações de alteração física.

```text
python -m pytest tests/test_diagnostic_capture.py tests/test_diagnostic_analysis.py tests/test_vorticity_budget.py tests/test_pressure_projection.py tests/test_lowmem_pressure.py -q
```

As verificações incluem passividade CPU/GPU, pressão de baixa memória, estiramento/inclinação analíticos com z não uniforme, telescopagem do orçamento, perda de identidade, RK4 em rotação exata e uma sequência sintética que percorre análise, trajetórias, orçamento e figuras.

Um teste adicional usa rotação sólida com w linear em z. A forma material contínua dá mudança Euleriana nula de zeta uniforme, enquanto o operador conservativo existente produz exatamente o termo adicional -zeta div(u) no interior. Esse caso analítico comprova por que o restante contínuo–discreto não pode ser chamado automaticamente de dissipação. Não demonstra que esse termo domina a execução da tempestade e não modifica o operador.
