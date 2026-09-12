# Pré-registro: influência causal do topo e do damping

Data do congelamento metodológico: 2026-09-11. Este documento antecede qualquer
integração completa do piloto. Resultados não serão usados para alterar critérios.

## Pergunta e escopo

O piloto pergunta se o topo em 15 km e seu damping impedem materialmente a
manutenção da convergência, do stretching e da concentração de vorticidade
vertical entre 0 e 2 km. Nenhum parâmetro de microfísica, LES, arrasto,
hodógrafa, termodinâmica, bolha ou intensidade do tornado será ajustado.

Este é um único par determinístico de triagem. Mesmo um resultado coerente não
será generalizado sem quatro pares independentes pré-registrados.

## Consolidação anterior

A proveniência v4 de 300 m terminou com gate `PASS`, 1.283 passos, neutralidade
bit a bit dos 12 campos prognósticos e máximo fechamento relativo
`8,7053e-15`. O HDF5 tem SHA-256
`fc8c5e3c57fad248294fa1e7c21328ae57d2fac1639eb320b14fcadfd10805b0`.

O par lateral 72/120 km rejeitou um grande efeito favorável na realização e
janela analisadas, mas não é uma prova causal geral. O caso de 120 km passou o
gate lateral; a generalização ainda exige ensemble.

## Auditoria do operador superior

`apply_velocity_bcs` é chamado uma vez na inicialização e quatro vezes por passo:
antes do preditor, depois do preditor, depois da projeção e depois do transporte
escalar/antes da microfísica. A chamada de inicialização não integra a contagem
por passo.

No controle de 48 níveis, o damping atua nas faces `45..48`, em
12.740,313, 13.457,105, 14.209,737 e 15.000 m. Os multiplicadores por chamada,
da face inferior à superior, são `0,95`, `0,977777...`, `0,994444...` e `1,0`.
A orientação é inversa à rampa absorvente usual: a maior redução fica na face
inferior e a face superior não é reduzida. O operador não usa `dt`; sua força
nominal por segundo depende do número de passos.

O projetor anelástico de baixa memória, usado nestas grades, fixa `w=0` nas faces
inferior e superior. Portanto, o tratamento efetivo é uma tampa rígida no teto
mais um sponge nas faces internas. A projeção pode repor `w` entre a segunda e a
terceira aplicação; a aplicação pós-projeção reintroduz divergência antes do
transporte escalar. Isso será medido, não reinterpretado como uma fronteira
aberta.

Serão separadas quatro quantidades:

1. força nominal definida pelos multiplicadores e número de chamadas;
2. alteração instantânea `Δw` e energia removida em cada chamada;
3. soma da perda realizada ao longo da trajetória integrada;
4. resposta causal do escoamento entre os dois braços.

O proxy offline anterior não será chamado de perda acumulada.

## Braços e única intervenção

- **CONTROL-15:** 72 × 72 × 15 km, `120×120×48`, configuração baseline;
- **HIGH-TOP-20:** 72 × 72 × 20 km, `120×120×54`, mesmas 49 faces entre
  0 e 15 km e seis novas células acima.

As novas faces são 15.829,776; 16.701,042; 17.615,870; 18.576,440;
19.585,039 e 20.000 m. Elas continuam geometricamente o espaçamento terminal;
a última célula é truncada para terminar exatamente em 20 km. As faces, centros
e volumes das primeiras 48 células permanecem bit a bit iguais ao controle.

O HIGH-TOP conserva quatro faces de damping e os mesmos multiplicadores. O
início da camada desloca-se para 17.615,870 m e sua espessura física é
2.384,130 m, contra 2.259,687 m no controle, diferença de 5,51%. Essa é a menor
diferença obtida ao continuar a malha e terminar exatamente em 20 km sem alterar
faces inferiores. A rampa invertida e as quatro aplicações por passo são
mantidas deliberadamente; corrigi-las exigiria executar novamente dois novos
controles e responderia outra pergunta.

O escalar legado `grid.dz`, usado como largura de filtro pela LES e por trechos
da microfísica, permanece `312,5 m` nos dois braços. Os operadores verticais,
sedimentação e volumes usam os `dz_c` locais. Alterar `grid.dz` para
`20000/54` mudaria o tratamento da coluna inferior e está proibido.

## Agenda temporal comum e abortamento

Os dois braços usarão os 5.421 passos nativos do controle arquivado em
`outputs/diagnostic_sequence_20260905/sequence.h5`, de `t=0` até
`3300,023494218 s`. A matriz little-endian `float64` das colunas
`t_start,dt,step_index` tem SHA-256
`6ebc5b14ac31b0f11819340e56331913cb233bd41453d617e4296501e742d5a3`;
`dt` varia entre 0,429081674 e 3,0 s.

Antes de cada passo, o runner calculará o limite adaptativo próprio do braço. Se
o `dt` imposto exceder esse limite por mais que tolerância de arredondamento, o
braço será interrompido antes do passo e a maior razão/violação observada será
documentada. Não será escolhida outra agenda silenciosamente. O preflight mede
o primeiro passo e verifica que `dz_min`, a largura de filtro e a coluna inferior
são idênticos; a segurança futura permanece um gate online porque a trajetória
HIGH-TOP ainda não existe antes do piloto.

## Instrumentação passiva

Um observador opcional, ausente por padrão, registra cada aplicação do damping:

- `w` antes por meio de estatísticas sem perda para o operador, `Δw` e `w` depois
  reconstruível exatamente como `w_antes+Δw`;
- energia removida por face com
  `0,5 ρ0_wface dx dy Δz_dual (w_antes²-w_depois²)`;
- perfis de `Δξ=∂Δw/∂y`, `Δη=-∂Δw/∂x` e `Δζ=0`;
- fluxo assinado `Σ p_dyn,w w dxdy` antes/depois;
- perfis de média e RMS de `w`, `p_dyn`, `p_dyn w` e convergência.

O volume dual da face superior é meia célula; faces internas usam metade da
soma das alturas adjacentes. Chamadas pré-projeção carregam a pressão dinâmica
anterior; chamadas pós-projeção carregam a pressão corrente. As quatro chamadas
do mesmo passo serão agregadas antes de correlações temporais.

Desde 2790 s, o observador também guarda, sem perda, `w_before` e `delta_w`
nas quatro faces de damping. O campo posterior é reconstruído exatamente pela
soma. Durante o spin-up, são gravados apenas perfis e integrais.

Reflexão não será declarada apenas por `p_dyn w`. Uma defasagem com progressão
vertical e reversão de fase é um proxy de propagação/reflexão; comprovação de
reflexão exigiria decomposição modal ascendente/descendente.

## Janela, rastreamento e métricas

Estados 3D completos e orçamentos serão capturados de 2790 a
3300,023494218 s nos mesmos 18 tempos-alvo e em volumes físicos comuns. Nenhum
campo será interpolado; apenas séries escalares poderão ser interpoladas dentro
do intervalo comum, sem extrapolação.

Métricas primárias:

1. `ζmax` a aproximadamente 121,7 m;
2. circulação assinada em raio de 4,2 km;
3. velocidade tangencial máxima;
4. largura de `ζ` à meia amplitude, com menor valor favorável;
5. convergência máxima de baixo nível;
6. stretching condicional médio de baixo nível.

Métricas secundárias: `ζmax`, integrais assinada e absoluta entre 0–2 km; RMW;
larguras de convergência e pressão; déficit de `p_dyn`; coerência e deslocamento
do eixo; tilting, dilatação, transporte, LES e projeção; energia do damping;
fluxos superiores; perfis de propagação e correlações com defasagem.

Cada métrica primária é materialmente favorável somente se tiver direção
favorável em pelo menos 12/18 tempos e mudança mediana de pelo menos 20%. Três
ou mais primárias materialmente favoráveis constituem resposta favorável
material do piloto. Uma ou nenhuma rejeita efeito favorável dominante neste
par. Dois resultados são classificados como mistos.

O rastreador deve permanecer no mesmo componente: saltos periódicos corrigidos
entre quadros devem ser menores que 6 km, e a correspondência entre braços deve
ser confirmada por continuidade da circulação, `w` e eixo vertical. Falha de
rastreamento torna a comparação `UNRESOLVED`.

## Gates congelados

HIGH-TOP-20 é geometricamente suficiente somente se todos forem satisfeitos:

- condensado total `>1e-5 kg/kg` não ocupa a última célula;
- nenhuma corrente ascendente `w>5 m/s` alcança a camada de damping;
- atividade `|w|>1 m/s` mantém margem física positiva abaixo do damping, com a
  margem mínima reportada em metros;
- todos os 12 campos prognósticos permanecem finitos;
- telescopagem do orçamento discreto de vorticidade tem erro relativo RMS
  máximo `≤1e-10`;
- observador ativo/inativo produz os 12 prognósticos bit a bit idênticos em CPU
  e GPU e registra exatamente quatro chamadas ordenadas por passo;
- o mesmo vórtice é rastreado de forma coerente.

Se a tempestade alcançar o damping, o teste será `GEOMETRICALLY INSUFFICIENT` e
não sustentará irrelevância do topo.

## Decisão e força da evidência

Se HIGH-TOP-20 passar todos os gates e não produzir resposta material, a
hipótese de que o topo é a causa dominante será enfraquecida para este caso. Se
produzir melhora material temporalmente associada a menor remoção ou menor
propagação descendente, o resultado sustentará desenhar um ensemble pareado de
quatro membros. O ensemble não será executado automaticamente.

Toda conclusão será classificada como: **comprovada pelos dados**, **fortemente
suportada**, **sugestiva** ou **não determinada**. O piloto termina antes de
qualquer mudança de física ou teste de resolução vertical perto da superfície.

## Correções após a primeira revisão independente

A primeira revisão somente leitura do Claude Code local não aprovou o início do
piloto até existirem: executor completo com abortamento CFL a cada passo, análise
e critérios congelados em código, e ensaio do fechamento v4 na grade de 54 níveis.
Esses três bloqueios foram corrigidos antes de qualquer integração completa.

Os dois braços usam os mesmos valores de referência `T_ref=266,84789304648075 K`
e `qv_ref=0,006483976752963724`. O preflight v3 mede o volume de saída a partir
das dimensões e da agenda, assume arrays sem compressão, acrescenta 25% para HDF5
e 256 MiB fixos, e exige ainda 3 GiB livres ao iniciar cada braço. O ensaio v4
de um passo na grade de 54 níveis obteve RMS relativo `3,29593e-16` para omega e
`3,62531e-16` para tilting. Os 15 testes específicos passaram. Uma segunda
revisão independente será registrada antes do uso da GPU.

## Produtos previstos

- `outputs/top_boundary_preflight_v3_20260911/summary.json`;
- `outputs/top_boundary_control15_20260911/sequence.h5` e
  `top_boundary_calls.h5`;
- `outputs/top_boundary_hightop20_20260911/sequence.h5` e
  `top_boundary_calls.h5`;
- `outputs/top_boundary_comparison_20260911/` com CSV, figuras, resumo e
  manifesto de hashes;
- `docs/TOP_BOUNDARY_CAUSAL_RESULTS.md` e atualização deste arquivo de
  continuidade.

