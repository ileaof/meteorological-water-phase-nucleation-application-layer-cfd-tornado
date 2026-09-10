# Experimento de extensão lateral da supercélula

Pré-registro: 2026-09-10, antes de existirem estados do período de análise.

O braço experimental altera somente `Lx=Ly` de 72 para 120 km e `nx=ny` de
120 para 200. Permanecem idênticos `dx=dy=600 m`, `Lz=15 km`, `nz=48`, todas as
faces verticais, sounding, movimento da tempestade, bolha térmica local, LES,
microfísica, arrasto, Coriolis, condições de contorno, damping e integração.
A janela analisada é 2790–3300,023494 s em tempos comuns, sem extrapolação.

O preflight exige que os únicos campos de configuração diferentes sejam as duas
extensões e as duas contagens horizontais, que os 12 campos prognósticos do
subdomínio central sejam bit a bit idênticos no início e que um passo real em
GPU permaneça finito.

As cinco métricas físicas primárias são zeta máxima a 121,7 m, circulação em
4,2 km, velocidade tangencial máxima, largura de zeta à meia altura e inclinação
do eixo entre a superfície e 2 km. Uma melhora material de uma métrica exige
efeito favorável em pelo menos 12 dos 18 tempos e mudança mediana de pelo menos
20% da magnitude mediana do controle. Três ou mais métricas assim classificadas
constituem efeito lateral favorável material; zero ou uma rejeitam a extensão
lateral como explicação favorável dominante; duas produzem classificação mista.

A suficiência geométrica requer que o condensado deixe de tocar a costura
periódica e que updrafts de 5 m/s permaneçam pelo menos 6 km da borda. Falhar
esses gates significa que 120 km ainda não capturou toda a tempestade, mesmo que
as métricas do vórtice mudem.

Este é um par determinístico de triagem, não um ensemble. Uma afirmação causal
forte exigirá os quatro pares planejados. O teste não avalia altura do topo nem
corrige a dependência do damping com o passo de tempo.
