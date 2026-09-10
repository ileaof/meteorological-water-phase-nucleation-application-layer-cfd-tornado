# Auditoria do tamanho do domínio — 9 de setembro de 2026

## Resultado

Há contato das estruturas condensadas com as bordas laterais periódicas e
penetração na região superior de amortecimento. O domínio de 72×72×15 km
não pode ser presumido independente das fronteiras. Isso não demonstra
que as fronteiras causam a deficiência de concentração de baixo nível.

Foram examinados os 18 estados completos de 2790,253 a 3300,023 s do caso
de 600 m, sem modificar nem reintegrar a solução. O teste de 300 m mantém
o mesmo domínio e não testa sensibilidade à extensão física.

## Medidas

Distâncias são dos centros de células até a borda geométrica. Uma distância
de 300 m significa ocupação da primeira/última célula horizontal.

| Indicador | Resultado na janela |
|---|---|
| Condensado total >0,01, >0,1 e >1 g/kg | Atinge células junto às bordas em todos os 18 estados |
| Área projetada com condensado >0,1 g/kg na faixa lateral de 6 km | 346,68–524,88 km² |
| Distância mínima das colunas com w>5 m/s às bordas | 0,3–7,5 km |
| Distância mínima das colunas com w>10 m/s às bordas | 5,1–9,9 km |
| Distância mínima do cold pool (anomalia de θv <−1 K, nível próximo de 100 m) | 11,1–12,9 km |
| Altitude máxima com condensado >0,01 g/kg | 13,83–14,60 km |
| Área com condensado >0,01 g/kg no último nível (14,60 km) | 0–50,4 km² |
| Altitude máxima com condensado >0,1 g/kg | 13,10–13,83 km |
| Máximo de w acima de 12 km | 5,92–9,28 m/s |
| Máximo de w nos centros acima da face inferior da região amortecida | 1,15–1,69 m/s |
| w na face superior nos estados salvos | zero |

Condensado significa a soma de ql, qi, qr, qs, qg e qh. Não é uma
identificação exclusiva de bigorna. Os limites incluem todas as estruturas,
inclusive possíveis células secundárias; não isolam apenas o vórtice rastreado.

## Condições de contorno e ressalva de implementação

Os metadados registram laterais periódicas e topo `damping_layer`.
Em `boundary_conditions.py`, o amortecimento atua nas últimas quatro faces
de w, começando em 12.740,31 m nesta grade. A implementação multiplica w
por fatores a cada aplicação de condição de contorno, sem escalar por dt;
o coeficiente é zero na face mais alta e cresce para baixo dentro da faixa.
Portanto, não deve ser descrita como uma camada Rayleigh convencional com
intensidade crescente até o teto e taxa temporal fixa. Há várias aplicações
por passo no núcleo. Esse detalhe merece auditoria própria e limita uma
interpretação do teste de 300 m como alteração exclusivamente espacial,
pois a mudança de passo temporal pode mudar o amortecimento efetivo.
Nenhuma alteração de física foi feita nesta auditoria.

## Interpretação e próximo teste

A suspeita sobre o domínio ganhou suporte geométrico, sobretudo nas laterais
para condensado e na vertical para a nuvem superior. O cold pool não toca
as bordas nesta janela. Contato com a costura periódica não prova interação
da tempestade consigo mesma nem identifica a causa da falha tornádica.

Uma comparação futura com 120×120×20 km deve preservar o espaçamento
horizontal e a grade vertical inferior, explicitar o tratamento do topo e
avaliar estruturas conectadas e fluxos nas fronteiras. Mudar domínio,
fronteiras e resolução simultaneamente não isolaria a causa. Não foi
iniciada nova simulação por esta auditoria.

## Reprodução

`python scripts/audit_domain_extent.py`

Fonte: `outputs/diagnostic_sequence_20260905/sequence.h5`.
Produtos: `outputs/domain_extent_audit/summary.json` e `timeseries.csv`.
