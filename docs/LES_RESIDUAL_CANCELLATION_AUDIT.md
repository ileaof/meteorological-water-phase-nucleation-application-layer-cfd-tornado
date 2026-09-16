# Auditoria do cancelamento no balanço LES

Data: 2026-09-15. Reanálise dos CSVs existentes, sem nova simulação.

## Conclusão

A conclusão anterior de que o transporte foi excluído como explicação do
inventário é forte demais. O cálculo publicado mede um resíduo líquido após
cancelamento temporal e o compara com produção integrada em módulo no espaço.
O resíduo permanece menor que essa produção, mas é comparável à produção
líquida dos intervalos. Isso não demonstra domínio do transporte nem uma causa
da ausência de tornadogênese; reabre uma atribuição antes considerada resolvida.

## Cálculo reproduzível

Em cada nível e intervalo, sejam `P` a produção LES assinada integrada no
disco, `G` a integral espacial do módulo do incremento LES e `R = delta I - P`
o resíduo arquivado na coluna `transport`.

- Medida anterior: `abs(sum R) / (sum G + abs(sum R))`.
- Sem cancelamento temporal de R: `sum abs(R) / (sum G + sum abs(R))`.
- Comparação com produção líquida por intervalo:
  `sum abs(R) / (sum abs(P) + sum abs(R))`.

Medianas entre níveis na janela 2790–3300 s:

| Medida | 600 m | 300 m |
|---|---:|---:|
| Medida anterior | 11,28% | 6,93% |
| Sem cancelamento temporal de R | 20,04% | 11,30% |
| Comparação com produção líquida por intervalo | 49,79% | 49,32% |
| Fração de R retida após soma assinada | 49,43% | 58,47% |

A 39,9 m, a terceira medida é 93,95% em 600 m e 97,87% em 300 m.
A 121,7 m, é 54,87% e 67,37%. São razões de magnitudes de termos do balanço,
não porcentagens causais da vorticidade final. A comparação com produção
líquida muda a pergunta e não substitui a comparação com G.

## Limitações

1. R não é uma medição independente de fluxo pela fronteira. No v4, o MUSCL
   transporta velocidades rotuladas e depois calcula seu rotacional. Interpretar
   toda essa evolução como entrada e saída de vorticidade exige decomposição adicional.
2. O módulo de R é tomado após integração espacial. Ainda há possibilidade de
   cancelamento entre células; G também contém cancelamento entre passos nativos.
3. A máscara é fixa nas duas pontas de cada diferença, mas muda entre intervalos.
   A soma não é a variação inicial-final em um único volume fixo.
4. O desfasamento de até 0,576 s em 600 m não fornece sozinho um piso de erro
   de 2%. É necessário medir ou limitar os incrementos nos trechos desencontrados.
   Os tempos em 300 m são coincidentes.
5. Mudança de sinal de um termo pequeno não demonstra ruído estatístico sem
   ensemble ou estimativa de erro. O resultado é sensibilidade entre dois casos.

Permanece sustentado que a atividade LES aumenta em toda a coluna e que o
maior aumento relativo não ocorre no primeiro nível. Não há nova justificativa
para a corrida de refinamento vertical cancelada.

## Teste de topo

Os metadados locais mostram CONTROL-15 completo: 5421 passos, 18 quadros,
tempo final 3300,023494 s. HIGH-TOP-20 foi abortado antes do passo de índice
1373, em 1287,331589 s, com zero quadros da janela de análise. Motivo:
`CFL_SCHEDULE_VIOLATION`, dt imposto 0,716317662 s contra limite próprio
0,716305057 s, razão 1,000017598. Não existe `sequence.h5` desse braço.

A comparação causal de topo permanece inconclusiva. Retomá-la exige revisar
explicitamente o protocolo de agenda comum e executar ambos os braços com
agenda admissível para ambos. Completar somente HIGH-TOP com outro dt não
satisfaz o isolamento original do damping por passo.

## Próximo diagnóstico

Separar a evolução do rótulo LES sob MUSCL de sua injeção local em intervalos
exatamente coincidentes e contabilizar o movimento da máscara. Para atribuir
entrada e saída, medir os termos de fronteira da formulação discreta escolhida.
O CSV atual não permite recuperar esses fluxos nem eliminar a defasagem de 600 m.

## Verificação e artefatos

Comando: `python scripts/audit_les_residual_cancellation.py`.

- Entrada: `outputs/les_local_vs_transported_20260913/les_local_vs_transported.csv`.
- Saídas: `outputs/les_residual_cancellation_20260915/by_height.csv` e `summary.json`.
- SHA-256 da entrada registrado no resumo.
- Identidade aritmética `delta I = P + R`: erro máximo zero no CSV.
- Caso sintético R=(+3,-3), P=(1,1), G=(2,2): razões 0, 0,6 e 0,75;
  verificação executada com sucesso.

A identidade verifica a reanálise, não valida independentemente a física.
