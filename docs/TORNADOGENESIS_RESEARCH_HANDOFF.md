# Tornadogênese: estado da pesquisa e passagem entre agentes

Atualizado em **2026-09-20**. Este documento é autocontido para leitura inicial.
Resume o estado deste projeto, não uma revisão geral da literatura científica.

## Resumo em um minuto

O modelo produz uma tempestade rotativa, mas os casos analisados não demonstram
um tornado resolvido. A descrição vigente é déficit de concentração/alinhamento
acompanhado de perda de circulação. Ainda não se identificou a causa física.

O avanço mais recente localizou uma contribuição negativa do estágio advectivo
MUSCL no balanço de circulação de baixo nível. Isso não demonstra dissipação
artificial: transporte e deformação físicos também passam por esse estágio.
A evolução do rótulo LES ainda é inferida por diferença, sem medição direta
dos seus fluxos. A próxima pergunta é como reconstruí-la a partir dos fluxos
de momento e distingui-la da injeção local e do movimento da região amostrada.

**Última etapa científica concluída:** pós-processamento em intervalos coincidentes.
**Etapa em preparação:** instrumentação passiva dos fluxos MUSCL.
**Ainda não executados:** piloto e replay diagnóstico novo de 600 m.

## Localização e escopo

- Projeto: `C:\Users\ileao\OneDrive\Documentos\met_h2o_nucleation_cfd_tornado`.
- Remoto: https://github.com/ileaof/meteorological-water-phase-nucleation-application-layer-cfd-tornado
- Branch de trabalho: `main`.
- Commit da análise discreta concluída: `f3460b6`.
- Janela madura: aproximadamente 2790–3300 s de tempo simulado.
- Malha de referência: 600 m horizontal; comparação existente: 300 m.
- Região diagnóstica: 17 níveis até aproximadamente 2 km; raios de 1,2; 2,4;
  4,2 e 6 km. O disco é centrado no rastreador de baixo nível.

O prompt anterior autorizou um replay de 600 m com piloto e gates. Ele não foi
iniciado. Este arquivo registra contexto e não é uma ordem de execução:
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
| evolução inferida do rótulo LES | +6.254 | +3.911 |
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

Pelo código v4, A é interpretado como evolução acumulada do rótulo sob MUSCL,
mas permanece inferido por diferença. Não é medição independente de fluxo
através da superfície. A soma por partes do rotacional já foi verificada;
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

**Não confundir:** a neutralidade anterior dos 12 prognósticos foi medida
entre replay observado e seu controle. A checagem nova contra o arquivo original
comparou zeta nas pontas, não os 12 prognósticos. O próximo replay deve fazer
explicitamente ambas as verificações. Fechamento aditivo não valida a física.

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
protege o comportamento existente, mas não valida o futuro observador completo.

Ainda faltam: observador completo conectado aos estágios; gravação das projeções
por passo; reconstrução independente pela divergência dos fluxos; executor
com agenda/gates; relatório do piloto; replay completo e sua interpretação.

## Próximos passos, na ordem

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

Referência histórica: 33 min para o replay de 600 m com controle. Estimativa
atual de armazenamento escalar: aproximadamente 65 MB brutos, ainda sem
medição pelo executor novo. Não anunciar tempo ou espaço como benchmarks.

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
| [Protocolo da captura](MUSCL_DIRECT_CAPTURE_PROTOCOL.md) | desenho da etapa ainda pendente |
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
