# Resultado da extensão lateral: 72 km versus 120 km

Data: 2026-09-10.

## Conclusão

O domínio de 120 x 120 x 15 km captura horizontalmente a supercélula durante
2790–3300 s pelo gate pré-registrado. A extensão lateral não melhora a
concentração do vórtice e pode ser rejeitada, neste par determinístico, como
explicação favorável dominante para a ausência de tornadogênese. O topo de
15 km continua alcançado pela nuvem e permanece uma limitação separada.

## Integridade e isolamento

O caso usa 200 x 200 x 48 células, `dx=dy=600 m`, a mesma grade vertical e a
mesma física do controle de 72 km. No preflight, somente `Lx`, `Ly`, `nx` e `ny`
diferiram. Os 12 campos prognósticos no subdomínio central foram bit a bit
idênticos no início. Um passo real em GPU passou com mais de 3,8 GB de VRAM
livre.

A integração terminou em 3300,023494 s após 5.578 passos. Foram validados 18
snapshots e 1.103 passos nativos na janela; todos os campos estavam finitos. O
snapshot final registra corretamente `step=5578`.

- HDF5: 4.993.784.422 bytes.
- SHA-256: `21169b14beb545d4e21698fe73b8d195e65b05272e71d4dfc4fa039584ca4513`.
- Tempo de parede: 4.433 s, aproximadamente 73,9 minutos.

## Métricas físicas comuns

Foram comparados 18 tempos comuns sem extrapolação, no nível de 121,7 m e em
volumes físicos iguais.

| Métrica primária | 72 km | 120 km | mudança pareada mediana | direção favorável |
|---|---:|---:|---:|---:|
| zeta máxima (s^-1) | 0,004739 | 0,004485 | -3,89% | 4/18 |
| circulação em 4,2 km (m2/s) | 20.157 | 19.146 | -5,89% | 5/18 |
| velocidade tangencial máxima (m/s) | 1,926 | 1,856 | -3,46% | 6/18 |
| largura de zeta à meia altura (m) | 2.142 | 2.294 | +4,53% | 6/18 |
| inclinação superfície–2 km (m) | 1.723 | 1.679 | -2,56% | 15/18 |

Nenhuma das cinco métricas alcançou simultaneamente o limiar de 20% e coerência
favorável em 12/18 tempos. O stretching condicional cresceu 10,9% em 18/18
tempos, mas permaneceu abaixo do limiar material. A zeta máxima entre 0 e 2 km
caiu 2,66%, a integral absoluta de zeta caiu 1,51% e a integral assinada caiu
4,57%. O sinal conjunto não indica que o domínio menor estivesse impedindo a
concentração.

## Bordas laterais e topo

| Distância mínima | 72 km | 120 km |
|---|---:|---:|
| condensado acima de 1e-5 | 0,3 km | 6,3 km |
| updraft acima de 5 m/s | 0,3 km | 26,7 km |
| cold pool abaixo de -1 K | 11,1 km | 35,7 km |

O caso de 120 km passa o gate lateral de 6 km, ainda que o condensado leve passe
por apenas 0,3 km de margem. No par e na janela analisados, 72 km era
geometricamente pequeno para a nuvem completa, mas a extensão não oferece
evidência de que essa proximidade explique o déficit do vórtice de baixo nível.

Condensado acima de `1e-5` ainda alcança o centro da última célula, a 14,605 km,
e ocupa até 41,04 km2 no último nível. Há movimento vertical de 1,05–1,34 m/s na
região de damping, iniciada em 12,740 km. Portanto, 15 km não permite declarar o
topo irrelevante.

## Damping e força da evidência

A trajetória de 120 km usa 1.103 passos na janela, contra 1.015 em 72 km. Como o
damping é aplicado quatro vezes por passo sem escala por `dt`, sua taxa nominal
na face inferior cresce de 0,4085 para 0,4439 s^-1, diferença de 8,65%. A
variável manipulada é a extensão lateral, mas a resposta total inclui esse
feedback numérico dependente do passo de tempo.

- **Comprovado neste par e janela:** em 120 km, `qcond>1e-5` ficou a pelo menos
  6 km da borda e as métricas de intensidade e circulação não aumentaram.
- **Fortemente suportado para esta realização:** a borda lateral de 72 km não
  era a causa dominante do déficit de concentração.
- **Não determinado:** efeito causal do topo/damping, suficiência de 20 km de
  altura e incerteza de ensemble.

Uma extensão vertical não deve alterar simplesmente `Lz` com `nz=48`, pois isso
engrossaria toda a grade inferior. O teste correto preserva as faces abaixo de
15 km e acrescenta níveis acima. Uma alegação causal forte sobre a extensão
lateral ainda exigiria os quatro pares planejados.

## Artefatos

- Preflight: `outputs/domain_extent_120km_preflight_20260910/summary.json`.
- Sequência: `outputs/domain_extent_120km_600m_20260910/sequence.h5`.
- Auditoria geométrica: `outputs/domain_extent_120km_audit_20260910/summary.json`.
- Comparação: `outputs/domain_extent_120km_comparison_20260910/summary.json` e
  `manifest.json`.
- Pré-registro: `docs/LATERAL_DOMAIN_EXPERIMENT.md`.
