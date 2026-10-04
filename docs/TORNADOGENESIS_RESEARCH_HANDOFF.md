# Tornadogênese: estado da pesquisa e passagem entre agentes

Atualizado em **2026-10-04**. Este documento é autocontido para leitura inicial.
Resume o estado deste projeto, não uma revisão geral da literatura científica.

## Resumo em um minuto

O modelo produz uma tempestade rotativa, mas os casos analisados não demonstram
um tornado resolvido. A descrição vigente é déficit de concentração/alinhamento
acompanhado de perda de circulação. Ainda não se identificou a causa física.

O avanço mais recente localizou uma contribuição negativa do estágio advectivo
MUSCL no balanço de circulação de baixo nível. Isso não demonstra dissipação
artificial: transporte e deformação físicos também passam por esse estágio.
A captura direta de 600 m confirmou a evolução do rótulo LES antes inferida
por diferença. Em 4,2 km e 121,7 m, x e y contribuem negativamente, enquanto
z compensa parte da perda; o MUSCL total é negativo nos 1015 passos nas três
convenções de máscara. São direções de fluxo de momento, não mecanismos
cinemáticos nem fluxos físicos de vorticidade demonstrados.

**Última etapa científica concluída:** captura direta MUSCL de 600 m e análise.
**Validação:** piloto e 1015 passos aprovados; 18 snapshots com os 12 campos
bit a bit iguais; observador neutro em todos os passos.
**Próxima etapa proposta, não executada:** decompor o fluxo nativo em referência
centrada e correção relativa em estados congelados, sem nova integração.
Relatório: [captura direta](MUSCL_DIRECT_CAPTURE_RESULTS.md).

**Revisão de 04/10:** auditoria independente dos arquivos aprovada, sem nova
simulação. No núcleo de 1,2 km/ponta final, operadores somam +4.378,74 m²/s,
mas máscara contribui -8.099,35, resultando em perda de inventário. A inversão
vertical de sinal depende da máscara, não define altura crítica física.
Detalhes e correções: [revisão dos dados](MUSCL_DATA_REVIEW_20261004.md).

## Localização e escopo

- Projeto: `C:\Users\ileao\OneDrive\Documentos\met_h2o_nucleation_cfd_tornado`.
- Remoto: https://github.com/ileaof/meteorological-water-phase-nucleation-application-layer-cfd-tornado
- Branch de trabalho: `main`.
- Commit da análise discreta concluída: `f3460b6`.
- Janela madura: aproximadamente 2790–3300 s de tempo simulado.
- Malha de referência: 600 m horizontal; comparação existente: 300 m.
- Região diagnóstica: 17 níveis até aproximadamente 2 km; raios de 1,2; 2,4;
  4,2 e 6 km. O disco é centrado no rastreador de baixo nível.

O replay de 600 m autorizado pelo prompt anexado foi concluído em 2026-09-20,
sem alterar física, malha ou agenda. Este arquivo não autoriza nova execução:
seguir o escopo da mensagem que ativar o próximo agente. Não alterar física,
amortecimento, LES, malha ou dt; não retomar o refinamento vertical cancelado.

## Vocabulário mínimo

| Termo | Significado no projeto |
|---|---|
| zeta | componente vertical da vorticidade, em s^-1 |
| circulação Gamma | integral horizontal assinada de zeta, em m²/s |
| LES local | incremento de velocidade efetivamente aplicado pelo operador subgrid |
| rótulo LES | parcela aditiva da velocidade, injetada pela LES e depois evoluída sob MUSCL |
| MUSCL | esquema de advecção de momento com reconstrução limitada nas faces da grade |
| máscara | conjunto de células incluídas no disco de análise; não é uma parcela material |
| fechamento | concordância entre a soma dos termos e a mudança observada |
| atribuição numérica | identificação de qual etapa do algoritmo alterou uma grandeza |
| causalidade física | identificação do mecanismo responsável por meio de controles apropriados |

As observações citadas abaixo são dos campos simulados. Não foi feita nova
validação meteorológica independente nesta etapa.

## O que se sabe

1. **Resolução:** 300 m intensifica e estreita o vórtice em relação a 600 m,
   mas ambos enfraquecem; o núcleo continua marginalmente resolvido. O par
   mistura malha, trajetórias e efeitos dependentes de dt. Não prova convergência.
2. **Distribuição radial:** a perda de circulação assinada é maior no anel
   de 1,2–4,2 km que no disco interno. Isso não prova exportação ou redução
   da vorticidade absoluta.
3. **Geometria:** em 300 m, o núcleo estreita e o eixo melhora entre as pontas
   enquanto a circulação cai. Enfraquecimento não equivale a espalhamento
   radial ou piora contínua do alinhamento.
4. **LES:** não foi demonstrada atuação preferencial junto ao solo nem causa
   LES para a ausência de tornado. O inventário rotulado não é produção local.
5. **Domínio:** ampliar as laterais para 120 km não trouxe melhora material
   neste par. HIGH-TOP-20 abortou em 1287,33 s, antes da janela madura; o
   teste causal do topo continua inconclusivo.
6. **História anterior:** o rótulo `initial` reúne tudo que precede o reinício.
   Não identifica sua origem física nem explica a gênese anterior a 2790 s.

## Resultado quantitativo principal

Acumulados em m²/s, a 121,7 m, disco móvel de 4,2 km, convenção de ponta final:

| Termo | 600 m | 300 m |
|---|---:|---:|
| MUSCL total | -19.312 | -22.516 |
| LES local | -5.667 | +869 |
| Coriolis | +7.103 | +6.087 |
| arrasto | -1.563 | +215 |
| movimento da máscara, circulação total | +162 | -13.513 |
| mudança da circulação total | **-19.277** | **-28.857** |
| evolução do rótulo LES (direta confirmada em 600 m; inferida em 300 m) | +6.254 | +3.911 |
| movimento da máscara, inventário LES | +2.756 | +257 |
| mudança do inventário LES | **+3.342** | **+5.037** |

Os demais operadores têm contribuição direta nula ou de arredondamento nessa
integral. Isso não elimina seus efeitos indiretos sobre a trajetória.

MUSCL total é negativo em 9/9 blocos nos dois casos em 4,2 km, tanto em
máscaras móveis (final e simétrica) quanto em disco fixo no centro inicial.
O sinal não é universal: em 1,2 km depende da convenção, e no disco fixo
de 6 km a soma MUSCL é positiva. Nove blocos não são nove realizações independentes.

No anel de 1,2–4,2 km, a perda é -15.556 e -27.704 m²/s; no disco interno,
-3.721 e -1.153 m²/s. O anel final de 300 m tem circulação assinada negativa;
razões entre circulações assinadas não são frações positivas limitadas a 1.

## Como o balanço foi definido

Seja C o rotacional discreto, q=C U_LES e I(M,q)=dx dy sum(M q).
Para duas pontas coincidentes a,b:

`L = C sum(delta U_LES local)`

`A = q_b - q_a - L`

`delta I = I(M_b,L) + I(M_b,A) + I(M_b-M_a,q_a)`

L e A são incrementos de vorticidade em s^-1. As três integrais têm m²/s.
O último termo mede mudança de seleção espacial. A convenção simétrica usa
a média das máscaras para L,A e a média dos campos para o termo de máscara.

Pelo código v4, A foi interpretado como evolução acumulada do rótulo sob MUSCL.
A captura de 20/09 confirmou essa interpretação em 600 m por incrementos
efetivos e reconstrução pelos fluxos; 300 m permanece sem essa captura direta.
Não é medição independente de fluxo físico através da superfície. A soma por
partes do rotacional já foi verificada;
ela representa a borda desse operador, não automaticamente fluxo físico de zeta.

Não somar essa atribuição por operadores com tilting/stretching cinemáticos:
são decomposições diferentes. Não chamar um resto de fechamento de difusão.

## O que foi efetivamente verificado

- Dez pontas coincidentes em 600 m delimitam nove blocos sem lacunas.
  Os mesmos índices foram agrupados em 300 m. Não se interpolaram campos.
- Zeta entre sequência e replay coincide exatamente nas pontas usadas.
- Soma dos operadores totais: erro RMS relativo máximo 8,07e-15.
- Identidades de máscara e soma por partes: erros absolutos inferiores
  a 5e-10 m²/s no diagnóstico concluído.
- Cinco testes analíticos do pós-processamento passaram.
- Manifestos registram hashes dos produtos, scripts e seleções HDF5 lidas.

**Atualização de 20/09:** a captura direta acrescentou comparação bit a bit dos
12 prognósticos em 18 snapshots originais, além da neutralidade por passo.
O maior erro relativo de fechamento por passo foi 5,58e-13, abaixo de 1e-10.
A diferença máxima da evolução LES direta para a inferida foi 5,46e-12 m²/s.
Fechamento aditivo não valida a física.

## Estado da implementação

As rotinas de momento e do tracer agora oferecem saídas opcionais para os
fluxos horizontais e verticais de momento u,v. A ordem das somas de tendência
foi preservada. `tests/test_momentum_flux_capture.py` verifica reconstrução e
passividade dessas rotinas isoladas, incluindo paredes, periodicidade e
ordens 1 e 2. Isso não substitui o piloto no estado maduro.

Em 20/09, os oito casos novos e os cinco testes do balanço discreto passaram
(13 testes). As verificações usam campos sintéticos em CPU, não o replay maduro.
Também passaram 29 testes existentes de proveniência e propagador, incluindo
passividade do tracer em CPU/GPU: **42 testes aprovados nesta etapa**. Isso
protege o comportamento existente, mas naquela etapa não validava o observador completo.

Depois desses testes, foram concluídos o observador, executor, piloto e replay.
Três testes adicionais do observador passaram, incluindo CPU/GPU. Os fluxos
são recalculados pelas rotinas nativas no estado anterior imutável; a mudança
efetivamente aplicada é verificada separadamente, com arredondamento explícito.
Não é uma validação independente da fórmula do esquema.

## Etapas concluídas em 20/09

1. **Fechar o protocolo.** Separar gates numéricos de hipóteses científicas;
   fixar escalas de tolerância por unidade; registrar se os fluxos foram
   capturados na chamada usada ou recalculados pela mesma rotina.
2. **Completar a captura.** Medir incrementos LES e MUSCL de total/rótulo;
   projetar contribuições x,y,z com pesos adjuntos; medir os outros estágios
   diretamente. Guardar módulos por passo e módulos de campos acumulados.
3. **Piloto curto.** Snapshot 93 até 94, 69 passos arquivados. Exigir agenda
   idêntica, neutralidade bit a bit e comparação dos 12 campos com o original;
   validar reconstrução, soma por partes e máscara. Medir custo real.
4. **Janela completa, após gates.** Mesmos 1015 passos até 3300,023494 s,
   verificando os snapshots disponíveis durante toda a execução. Falha numérica
   interrompe a integração; sinal inesperado de um termo não é falha numérica.
5. **Análise por direção e escala.** Comparar captura direta com A inferido;
   examinar discos e anéis, máscaras e cadências. Esclarecer o que é robusto
   antes de propor qualquer intervenção física ou estudo de convergência.

Todas as cinco etapas acima foram executadas. Custo medido do novo replay:
21,43 min e 44,27 MiB de HDF5, sem novos campos completos. Preservar a lista
como registro do que foi feito, não como ordem para repetir a execução.

## Próximo passo mínimo

Pós-processar os 18 estados arquivados, sem integração: definir um fluxo
centrado de referência nas mesmas faces e a correção MUSCL menos referência;
verificar a reconstrução e comparar suas projeções assinadas e módulos nos
mesmos discos/anéis. Essa diferença inclui reconstrução e upwinding: não isola
o limitador nem identifica automaticamente erro ou dissipação artificial.
Uma referência upwind de primeira ordem pode ajudar a distinguir essas parcelas.
Estados arquivados não são todos os estados pré-MUSCL: não tratar a amostragem
esparsa como integral exata da janela. Protocolo, critérios de interpretação
e estimativa de custo estão no relatório novo. Nenhum novo replay está autorizado
por este documento.

## Resultado novo em números

Em 121,7 m, disco móvel de 4,2 km, ponta final, m²/s: x=-24.769,69;
y=-22.388,27; z=+27.846,31; MUSCL=-19.311,65. O total é negativo nos
1015 passos nas três convenções. No anel 1,2--4,2 km, MUSCL=-21.785,77;
em 4,2--6 km, +10.168,67. O sinal depende de raio/altura: não há perda
universal. O cancelamento espacial do MUSCL total é forte, mas não há
cancelamento temporal da integral assinada nesse disco de referência.

## Perguntas que essa captura pode responder

- A soma dos fluxos medidos reproduz o incremento MUSCL e a evolução do rótulo?
- Quais direções de fluxo de momento e regiões contribuem mais para a perda?
- Quanto se cancela dentro de cada passo, entre passos e entre regiões?
- Quais conclusões mudam com raio, máscara e cadência?

Ela não resolve sozinha a separação entre dinâmica física e erro de discretização.
Mesmo uma direção vertical dominante não demonstra tilting, estiramento ou
exportação vertical de zeta sem uma derivação discreta adicional.

## Conclusões substituídas: não reutilizar

- "O inventário LES é demonstradamente produzido localmente e transporte foi excluído."
- "Maior quota LES demonstra ausência de convergência ou atuação preferencial superficial."
- "O par 300/600 parte do mesmo estado maduro e só dx muda."
- "Mudança de sinal é ruído estatístico", sem erro estimado ou ensemble.
- "Desfasamento de 2% implica erro de 2%", sem medir incrementos desencontrados.
- "Fechamento próximo à precisão de máquina valida o modelo físico."
- "O teste de topo terminou", deduzido apenas da existência de sua pasta.

Textos antigos foram preservados como histórico; priorizar as atualizações de
15/09 e este documento de 20/09. Não usar instruções antigas como ordens atuais.

## Mapa de leitura e reprodução

| Arquivo | Uso |
|---|---|
| [Resultados discretos](LES_DISCRETE_BALANCE_RESULTS.md) | relatório científico completo, tabelas e limitações |
| [Método discreto](LES_DISCRETE_BALANCE_METHOD.md) | definições matemáticas e unidades |
| [Auditoria de cancelamento](LES_RESIDUAL_CANCELLATION_AUDIT.md) | correção da interpretação de resíduos |
| [Protocolo da captura](MUSCL_DIRECT_CAPTURE_PROTOCOL.md) | desenho e critérios da etapa concluída |
| [Resultados da captura](MUSCL_DIRECT_CAPTURE_RESULTS.md) | direções, máscaras, cancelamentos e próximo diagnóstico |
| [Continuidade](TORNADOGENESIS_ANALYSIS_CONTINUITY.md) | histórico datado; contém propostas antigas |
| [Análise por altura](CLOSURE_CONVERGENCE_BY_HEIGHT.md) | contraste LES e correções sucessivas |

Scripts concluídos: `scripts/analyze_les_discrete_balance.py`,
`scripts/summarize_les_discrete_balance.py`, `scripts/audit_les_residual_cancellation.py`.
Eles fazem pós-processamento; seus comandos não são lançadores de simulação.

Resultados locais: `outputs/les_discrete_balance_20260915_v2/` e
`outputs/les_discrete_balance_synthesis_20260915/`, com manifestos.
Cópia compacta versionada: `docs/media/storm/les_discrete_balance_20260915/`.
Os grandes HDF5 são locais e não estão incluídos no remoto:

- `outputs/diagnostic_sequence_20260905/sequence.h5`;
- `outputs/vorticity_provenance_long_v4_2790_3300/provenance_long_v4.h5`;
- `outputs/resolution_300m_20260909/sequence.h5`;
- `outputs/resolution_300m_provenance_20260909/provenance_long_v4.h5`.

Um agente com apenas o remoto pode revisar conclusões e tabelas compactas,
mas não reproduzir o processamento dos campos nem executar o replay sem esses
arquivos. Não substituir arquivos ausentes por uma simulação desde t=0.

## Regra para a próxima passagem

Registrar separadamente: pergunta examinada, dados usados, mudanças de código,
testes executados, simulações realmente iniciadas, resultados, limitações e
próxima decisão. Acrescentar uma entrada datada à continuidade somente com
resultados ou mudanças de estado verificáveis. Nunca relatar uma proposta
como execução concluída.
