# Auditoria de proveniência da vorticidade horizontal

## Resposta executiva

**Os dados atuais não permitem reconstruir com segurança a proveniência de `ωh` sem nova integração.** Eles permitem calcular `ξ`, `η`, seus balanços intervalares e o curl de cada incremento de velocidade, mas os incrementos estão somados em janelas de aproximadamente 30 s. A identidade de origem precisa ser transportada e deformada em cada passo nativo, na ordem efetiva dos operadores. Essa propagação não pode ser invertida a partir das somas arquivadas.

Foi implementada a instrumentação passiva mínima para preencher essa lacuna. Ela:

- mantém nove contribuições vetoriais de vorticidade no domínio completo;
- aplica a todas elas o mesmo transporte, stretching, reorientação e dilatação resolvidos;
- injeta em cada origem o curl exato do incremento nativo de seu operador;
- guarda separadamente a diferença entre a advecção nativa em forma de fluxo e a reconstrução cinemática contínua;
- salva somente `ξ_j`, `η_j`, `ξ`, `η` e `∇h w` entre 0 e 2 km mais três níveis de halo;
- não escreve em qualquer variável prognóstica.

A instrumentação passou em testes CPU e GPU, foi bit a bit neutra em execuções pareadas e fechou `ω` e `Tz` próximo do arredondamento. Esses testes validam a instrumentação; não fornecem ainda uma atribuição física para a tempestade.

**Classificação final: REQUER NOVA INSTRUMENTAÇÃO PASSIVA.** A instrumentação necessária já está implementada e validada em execução curta. A execução longa não foi iniciada.

## 1. Equação efetivamente diagnosticada

O solver prognostica `u`, `v` e `w` na malha C de Arakawa. A vorticidade é uma variável diagnóstica em centros de célula:

`ξ = ∂w/∂y − ∂v/∂z`,

`η = ∂u/∂z − ∂w/∂x`,

`ζ = ∂v/∂x − ∂u/∂y`.

Para a velocidade resolvida e para cada força/incremento `F`, a identidade contínua correspondente é:

`Dω/Dt = (ω·∇)u − ω(∇·u) + curl(F)`.

Em componentes:

`Dξ/Dt = ξ ∂u/∂x + η ∂u/∂y + ζ ∂u/∂z − ξ ∇·u + [curl(F)]x`,

`Dη/Dt = ξ ∂v/∂x + η ∂v/∂y + ζ ∂v/∂z − η ∇·u + [curl(F)]y`,

`Dζ/Dt = ξ ∂w/∂x + η ∂w/∂y + ζ ∂w/∂z − ζ ∇·u + [curl(F)]z`.

Na terceira equação, `ξ∂w/∂x + η∂w/∂y` é o tilting vertical, `ζ∂w/∂z` é o stretching axial e `−ζ∇·u` é a dilatação. Como o modelo é anelástico, `∇·u` não é imposto igual a zero; a restrição é `∇·(ρ0u)≈0`. Por isso o termo de dilatação deve ser mantido.

Esta forma contínua não substitui o operador real. A advecção de momento implementada é:

`∂u_i/∂t = −∂(u_j u_i)/∂x_j`,

em forma de fluxo, com reconstrução MUSCL/minmod de segunda ordem nos volumes escalonados. O diagnóstico preserva explicitamente a diferença entre o curl desse incremento nativo e a identidade contínua discretizada.

## 2. Mapeamento para o solver real

A ordem executada em cada passo é: contornos iniciais, LES, advecção de momento, buoyancy, Coriolis, drag, forçantes externas, guard, contornos do preditor, projeção, contornos da projeção, transporte escalar, microfísica e contornos finais.

| Estágio real | Implementação | Ação sobre `ωh` | Classe |
|---|---|---|---|
| contornos iniciais | `core._predictor`, `apply_velocity_bcs` | substituição/cópia de velocidades nas faces; curl do incremento é retido como `boundary` | A/B, dependente do contorno |
| LES | `strain_and_viscosity` e `apply_les_momentum` | `div(Km grad u_i')`; o curl nativo pode criar ou destruir contribuição assinada | A e efeito dissipativo |
| advecção | `add_momentum_advection` | transporte de `ω`, stretching, rotação entre componentes e dilatação; diferença discreta explícita | B/C/D |
| buoyancy | `buoyancy_w_tendency`; `w += dt Bf` | força `(0,0,B)` interpolada nas faces verticais | A direta em `ωh` |
| Coriolis | `du/dt=fv'`, `dv/dt=−fu'` | curl do incremento horizontal; pode gerar e reorientar vorticidade | A/D |
| surface drag | lei implícita de tensão superficial | curl da desaceleração nos níveis dentro da camada superficial | A e E |
| external forcing | aquecimento/umidificação opt-in | no código atual muda escalares, não momento; atua depois via `B` | E |
| guard | clip extremo de `u,v,w` | se ativado, é uma alteração numérica direta; armazenada em `other` | A numérica |
| predictor/projection/final BCs | `apply_velocity_bcs` | alteração direta de faces, armazenada em `boundary` | A/B |
| projection | correção anelástica `−dt (1/ρ0_face) grad(p')` | curl do incremento exato; o termo vertical contínuo tende a zero, mas o horizontal não precisa ser zero quando `ρ0=ρ0(z)` | A/D e E |
| transporte escalar e microfísica | advecção/difusão de escalares, conversões e sedimentação | não alteram diretamente a velocidade nesse estágio; mudam a buoyancy dos passos seguintes | E |

Aqui, A significa geração/destruição direta, B transporte, C amplificação/stretching, D reorientação e E efeito indireto. LES, drag e projection são armazenados como contribuições assinadas do processo. Um valor negativo representa remoção/cancelamento de vorticidade criada por outras origens; não significa uma nova espécie material de vorticidade.

Para a buoyancy efetivamente usada,

`B = g[θ'/θ0 + 0,61 qv' − (ql+qi+qr+qs+qg+qh)]`,

quando `moisture_buoyancy` está ativa. A força é aplicada somente a `w`, portanto:

`curl(B k) = (∂B/∂y, −∂B/∂x, 0)`.

Isso é geração de vorticidade horizontal pelo curl da força de flutuabilidade implementada. Não é identificado aqui como a forma baroclínica compressível completa.

## 3. Auditoria dos HDF5 existentes

| Arquivo | Quadros | Passos nativos | Conteúdo relevante |
|---|---:|---:|---|
| sequência histórica 0–3900 s | 131 | 6414 | estado 3D, incrementos por operador até 2 km + halo e integrais cinemáticas totais |
| cada ramo causal 2370–3300 s | 32 | 1985 | estado e incrementos somente nos 20 níveis inferiores, mais integrais cinemáticas totais |

Nas primeiras janelas históricas há 11–13 passos nativos por quadro. Nos ramos causais maduros há aproximadamente 68 passos por quadro. Cada grupo `increments/<operador>` contém a soma das mudanças de `u`, `v` e `w` durante toda a janela. O grupo `kinematic_integrals` contém integrais da vorticidade total, não a evolução de contribuições rotuladas.

Faltam exatamente:

1. `ξ_j`, `η_j` e `ζ_j` no início de cada passo;
2. o gradiente completo de velocidade em cada passo, na ordem em que atua sobre cada origem;
3. o curl de cada incremento antes de ser agregado com dezenas de passos posteriores;
4. a composição por origem da vorticidade importada verticalmente de níveis acima de 2 km;
5. a proveniência anterior ao reinício causal em 2370 s.

Somar curls históricos numa célula não resolve o problema: as contribuições criadas no início de uma janela são transportadas e giradas antes do fim, e os operadores não comutam. Também não é possível aplicar retrospectivamente um único tensor de deformação médio a todas as fontes, pois cada uma nasce em tempo e posição diferentes.

Consequentemente, os HDF5 podem sustentar um orçamento instantâneo/intervalar de `dξ/dt` e `dη/dt`, mas não a pergunta “onde essa vorticidade nasceu e no que se transformou depois?”.

## 4. Formulação dos traçadores passivos

Para cada rótulo `j`, a instrumentação evolui:

`∂ω_j/∂t = −u·∇ω_j + (ω_j·∇)u − ω_j∇·u + S_j`.

O operador homogêneo é linear em `ω_j` quando a velocidade resolvida é fixada. Assim, todas as origens experimentam o mesmo transporte, stretching, reorientação e dilatação. As fontes são adicionadas na ordem real do split do solver:

`Δω_j,source = curl(u_after,j − u_before,j)`.

Para a advecção:

`Δω_adv,native = curl(Δu_adv,native)`,

`R_adv = Δω_adv,native − dt Σ_j L_u(ω_j)`.

`R_adv` é acumulado sob `advection_remainder`. Ele inclui a diferença entre forma de fluxo e forma contínua, limitador, interpolação C-grid e discretização. Não é chamado automaticamente de difusão numérica.

Os rótulos implementados são:

- `initial`: vorticidade existente quando o traçador é ativado;
- `buoyancy`;
- `les`;
- `surface_drag`;
- `coriolis`;
- `projection`;
- `boundary`;
- `other`, incluindo external forcing direto e guard;
- `advection_remainder`.

Para uma execução definitiva iniciada em `t=0`, `initial` representa a vorticidade do hodógrafo/estado inicial. Se a instrumentação for ativada num restart, o nome deve ser interpretado como **antecedente no instante do restart**; sua origem anterior permanece indeterminada.

## 5. Tilting e geometria por origem

Nos quadros sincronizados serão reconstruídos:

`Tz,j = ξ_j ∂w/∂x + η_j ∂w/∂y`,

`Tz,j = |ωh,j| |∇h w| cos(φ_j)`,

`φ_j = atan2/acos` consistente entre `ωh,j` e `∇h w`.

O fechamento é:

`R_T = Tz,total − Σ_j Tz,j`.

Como o produto é linear em `ωh,j` para o mesmo `∇h w`, `R_T` é determinado diretamente por `Rξ` e `Rη`. Serão tratados separadamente magnitude, orientação e contribuição assinada. “Fração por origem” não será calculada ingenuamente como uma porcentagem positiva quando fontes se cancelarem; serão mostradas magnitudes, projeções assinadas e a soma algébrica.

## 6. Implementação mínima realizada

A implementação está em `src/storm_dynamics/vorticity_provenance.py`. Ela usa os hooks `diagnostic_observer` já existentes, sem alterar `StormSimulation`, o integrador ou qualquer configuração física.

Características:

- campos internos `ω_j=(ξ_j,η_j,ζ_j)` no domínio vertical completo;
- operações no mesmo backend NumPy/CuPy do solver;
- uma origem processada por vez para limitar temporários;
- HDF5 separado, em precisão `float64` por padrão;
- saída de `ωh,j`, `ωh,total` e `∇h w` apenas entre 0–2 km mais três halos;
- histórico de fechamento em todos os estágios;
- checkpoint completo dos três componentes para reinício/fork;
- `DiagnosticObserverChain` para uso simultâneo com a captura 3D já existente.

Manter os três componentes no domínio completo durante a integração é necessário. `ξ_j` e `η_j` recebem contribuição de `ζ_j∂u/∂z` e `ζ_j∂v/∂z`, e ar acima de 2 km pode voltar à camada analisada. Cortar os traçadores em 2 km faria a proveniência da importação vertical virar uma condição de contorno desconhecida. O corte é aplicado somente na escrita.

## 7. Neutralidade e fechamento

Foram executadas duas validações curtas.

### Teste pequeno pareado

- malha `10×10×8`, seis passos de 0,1 s;
- todos os 16 campos comparados (`u,v,w,p,θ,qv`, seis hidrometeoros e diagnósticos termodinâmicos) foram bit a bit idênticos ON/OFF;
- erro máximo em todos os prognósticos: exatamente zero;
- `max|Rω| = 5,20×10⁻¹⁸ s⁻¹`;
- `RMS(Rω)/RMS(ω) = 1,52×10⁻¹⁶`;
- `max|RT| = 1,06×10⁻²² s⁻²`;
- `RMS(RT)/RMS(Tz) = 1,98×10⁻¹⁶`;
- checkpoint/restart das fontes: bit a bit idêntico.

### Teste no estado maduro

- reinício no quadro 93 da sequência histórica, `t=2790,253490 s`;
- malha real `120×120×48`, GPU, dois passos nativos arquivados de 0,43546 e 0,43551 s;
- os 12 campos prognósticos arquivados foram bit a bit idênticos ON/OFF;
- erro máximo prognóstico: zero;
- `max|Rω| = 1,39×10⁻¹⁷ s⁻¹`;
- `RMS(Rω)/RMS(ω) = 1,19×10⁻¹⁶`;
- máximo de `|RT|` entre todos os estágios: `8,13×10⁻²⁰ s⁻²`;
- `RMS(RT)/RMS(Tz) ≤ 1,63×10⁻¹⁶`.

Os testes automatizados específicos passaram em CPU e GPU. A bateria combinada de proveniência e diagnósticos existentes passou com 16 testes.

Para a sequência longa, a atribuição somente será aceita se, simultaneamente:

- `RMS(Rωh)/RMS(ωh) ≤ 10⁻¹⁰`;
- `max(|Rξ|,|Rη|) ≤ 10⁻¹⁰ s⁻¹`;
- `RMS(RT)/RMS(Tz) ≤ 10⁻¹⁰`;
- `max|RT| ≤ 10⁻¹² s⁻²`;
- o residual for menor que 1% da menor contribuição usada numa conclusão física.

Os limites são muito mais frouxos que o arredondamento medido, mas muito menores que as escalas de `ωh`, `Tz` e dos contrastes anteriores. Se qualquer limite falhar ou se `advection_remainder` for comparável às origens interpretadas, a atribuição será interrompida.

## 8. Custo de memória, VRAM e disco

Para `120×120×48`, nove origens, três componentes internos e `float64`:

| Item | Custo |
|---|---:|
| campos persistentes do instrumento, incluindo `ω_total` e cópia da velocidade | 182.684.160 bytes = 174,2 MiB |
| campos de origem sem auxiliares, para checkpoint completo | 142,4 MiB |
| quadro de saída 0–2 km + halo, sem compressão | 50.688.000 bytes = 48,34 MiB |
| 32 quadros de um ramo | 1,51 GiB brutos |
| 32 quadros em três ramos | 4,53 GiB brutos |

Em GPU, os 174,2 MiB persistentes residem em VRAM. Os temporários de derivadas são reutilizados origem a origem; recomenda-se reservar aproximadamente 350–400 MiB adicionais de VRAM para o instrumento durante o passo, a confirmar por benchmark do alocador antes da execução longa. Na escrita, um dataset por vez é transferido ao host, evitando uma cópia simultânea de todo o conjunto.

No teste maduro, dois passos custaram 0,699 s sem instrumento e 0,726 s com instrumento, razão 1,039. Essa amostra é curta e inclui sincronizações de métricas; serve como indicação, não como benchmark definitivo. No teste minúsculo CPU, a razão foi 2,36 porque o overhead fixo domina. Antes da execução longa deve-se medir 30–60 s do estado maduro com e sem escrita.

O HDF5 usa gzip nível 1 e shuffle. O tamanho comprimido depende dos campos e não é usado para garantir espaço. O planejamento deve reservar o tamanho bruto, além do arquivo diagnóstico convencional.

## 9. O que precisa ser salvo

Não é necessário salvar dezenas de campos 3D completos em cada passo. O conjunto mínimo por quadro é:

- `ξ_j`, `η_j` para as nove origens, 0–2 km + halo;
- `ξ_total`, `η_total` no mesmo recorte;
- `∂w/∂x`, `∂w/∂y` no mesmo recorte;
- métricas de fechamento e metadados.

`Tz,j`, ângulos, magnitudes e orientações são reconstruídos offline desses campos. As trajetórias e máscaras Eulerianas usam a sequência 3D sincronizada convencional. Internamente, os três componentes `ξ_j,η_j,ζ_j` permanecem no domínio completo. Um checkpoint de fork precisa salvar todos eles uma vez.

## 10. Execução curta mínima

A validação mínima já executada cobre:

1. todos os estágios do split por seis passos;
2. identidade bit a bit ON/OFF em CPU;
3. o mesmo teste automatizado em GPU;
4. fechamento de `ω` e `Tz` em todos os estágios;
5. escrita e leitura bit a bit do checkpoint;
6. dois passos na malha real e num estado maduro próximo da inversão de tilting.

Antes de autorizar a sequência longa, resta somente um benchmark de 30–60 s no estado maduro se for necessário confirmar tempo total e pico real do alocador. Ele não é necessário para validar a matemática ou a neutralidade já demonstradas.

## 11. Execução longa necessária

A execução cientificamente completa deve preservar a origem anterior a 2370 s:

1. integrar uma vez o caso compartilhado de `t=0` a aproximadamente 2370 s com os traçadores ativos, sem salvar quadros de proveniência frequentes;
2. gravar um checkpoint completo de estado prognóstico e de `ω_j` em 2370 s;
3. clonar esse checkpoint para WEAK-EVAP, CONTROL e STRONG-EVAP;
4. aplicar somente os fatores já definidos no experimento causal, a partir do mesmo instante de intervenção;
5. integrar os três ramos de 2370 a 3300 s, salvando a cada 30 s os recortes de proveniência e a sequência diagnóstica convencional;
6. analisar detalhadamente 2790, 2820, 2850, 2910, 2940 e 3000 s e as 27 trajetórias equivalentes.

Ativar o traçador apenas no restart de 2370 s é uma alternativa menor e responde quanto da diferença pós-intervenção vem de novas fontes. Ela deixa toda a vorticidade anterior agrupada em `initial/antecedent` e não responde plenamente de onde veio a `ωh` presente na inversão. Por isso não é a execução definitiva recomendada.

Não é necessário integrar até 3900 s para a pergunta atual; 3300 s cobre a janela solicitada. Nenhuma dessas integrações longas foi executada nesta etapa.

## 12. Análise futura após a sequência

Com os novos dados será possível produzir, para WEAK, CONTROL e STRONG:

- `|ωh,j|`, orientação, ângulo com `∇h w` e `Tz,j` nas 27 trajetórias;
- contrastes SW, SC e CW de magnitude, posição, orientação e tilting por origem;
- médias, máximos, integrais e mapas por origem no componente Euleriano;
- distâncias entre cada origem, `ζmax`, circulação, `Tz`, `|∇h w|`, cold pool e `|∇h B|`;
- resíduos `Rξ`, `Rη` e `RT` no tempo e no espaço;
- as 18 figuras de proveniência solicitadas.

Essas figuras físicas não foram produzidas agora porque os dados que as sustentariam ainda não existem. Foi produzida apenas a figura de fechamento da validação curta.

## 13. Respostas às dez perguntas

1. **Os dados atuais permitem reconstruir a proveniência sem nova simulação?** Não. Os operadores foram agregados em janelas de 30 s e não há estado de origem por passo.
2. **Quais informações faltam?** Campos `ω_j`, gradiente de velocidade e curl dos incrementos na ordem de cada passo, inclusive acima de 2 km e antes do restart causal.
3. **Qual é a formulação correta?** Uma decomposição vetorial que evolui cada `ω_j` com o operador cinemático linear comum e injeta fontes por seus curls nativos, com resto advectivo explícito.
4. **Quais campos foram adicionados?** Nove vetores internos completos; na saída, somente os componentes horizontais, `ωh,total`, `∇h w` e métricas.
5. **Como garantir passividade?** O módulo usa somente cópias e leituras pelo observer; nenhuma referência retornada alimenta o solver. A identidade ON/OFF foi testada bit a bit.
6. **Como verificar os fechamentos?** Calcular `Rξ`, `Rη` depois de todos os estágios e `RT` usando o mesmo `∇h w`; aplicar as tolerâncias quantitativas acima.
7. **Qual o custo?** 174,2 MiB persistentes, temporários adicionais estimados em até 350–400 MiB e 48,34 MiB brutos por quadro de proveniência.
8. **São necessários campos 3D completos?** Durante a evolução, sim. Na saída periódica, não: bastam os componentes horizontais abaixo de 2 km mais halo. Um checkpoint completo é necessário no fork.
9. **Qual execução curta valida?** Seis passos pequenos CPU/GPU mais dois passos GPU na malha real e no estado de 2790 s; ela já foi concluída e passou.
10. **Qual execução longa será necessária?** Um trecho compartilhado 0–2370 s, checkpoint/fork e três ramos 2370–3300 s. Ela aguarda autorização.

## 14. Artefatos

- `src/storm_dynamics/vorticity_provenance.py`: traçadores, fechamento, HDF5 e checkpoint;
- `tests/test_vorticity_provenance.py`: linearidade, fonte de buoyancy, passividade CPU/GPU, fechamento e restart;
- `scripts/validate_vorticity_provenance.py`: validação curta pareada;
- `scripts/validate_vorticity_provenance_mature.py`: validação de dois passos na malha real;
- `outputs/vorticity_provenance_validation_v4/summary.json`: resultados da validação curta;
- `outputs/vorticity_provenance_validation_v4/closure_residuals.png`: resíduos por estágio;
- `outputs/vorticity_provenance_validation_mature/summary.json`: resultados do estado maduro.

## Limites de interpretação

Os rótulos são uma decomposição aditiva de processo e origem, não experimentos contrafactuais. Mesmo que `Tz,buoyancy` venha a dominar, a conclusão permitida será que a vorticidade atribuída ao curl da força de buoyancy participa do tilting. Isso, isoladamente, não prova que buoyancy causou tornadogênese.

O rótulo `advection_remainder` deve permanecer visível. Se ele ou os resíduos de fechamento forem comparáveis a uma origem candidata, não haverá atribuição robusta para essa origem. Boundary e projection também devem ser interpretados no contexto do operador anelástico e da malha escalonada.

## Classificação final

**REQUER NOVA INSTRUMENTAÇÃO PASSIVA.**

A instrumentação passiva mínima está implementada, é reiniciável, fecha numericamente e foi bit a bit neutra. Os arquivos históricos não contêm a história de deformação por origem necessária. A próxima ação científica seria a execução longa descrita acima, somente após autorização explícita.
