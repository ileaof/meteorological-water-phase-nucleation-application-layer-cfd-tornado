# Precursor de tornado: da regra verbal a uma quantidade numérica

Data: 2026-09-13.

## A pergunta

«Vento à superfície num sentido, a nuvem a cisalhar no sentido oposto, e humidade
a subir a partir do solo — como transformo isso numa previsão numérica?»

## A tradução direta

A descrição é, em rigor, uma afirmação sobre **helicidade**: um hodógrafo que
gira. As três peças numéricas são:

1. **Vorticidade horizontal ambiente** — `ω_h = (−∂v/∂z, ∂u/∂z)`. «Sentidos
   opostos» entre o solo e a camada da nuvem é |ω_h| grande.
2. **Fração streamwise** — a projeção de ω_h no vento *relativo à tempestade*,
   `ω_s = ω_h·(V−C)/|V−C|`. Só a parte streamwise gera um updraft **rotativo**;
   a crosswise produz um par que se divide, não uma rotação. Note-se que ω_h é
   sempre perpendicular ao vetor de cisalhamento, pelo que vorticidade
   streamwise exige o escoamento relativo à tempestade **perpendicular** ao
   cisalhamento — é isso que a curvatura do hodógrafo fornece.
3. **SRH** — `∫ (v−c_y)∂u/∂z − (u−c_x)∂v/∂z dz`, em 0–1 km e 0–3 km, com C de
   Bunkers.

A humidade a subir entra por outra via: **LCL baixo**, CAPE, CIN.

## A correção que as medições deste projeto obrigam

Estes três ingredientes são **necessários e não suficientes**, e aqui isso não é
uma cautela retórica — foi medido:

- **Tentativa E**: subir a SRH de 254 para 648 m² s⁻² mudou a rotação de baixo
  nível em ~0%. A SRH foi eliminada como limitante.
- **Tentativa H**: `vorticity_budget.tilting_efficiency` fatoriza
  `inclinação = |ω_h| · |∇_h w| · cos θ`. Nos campos do modelo ω_h era abundante
  (1,1e−2 s⁻¹ a 51 m) e ∇_h w existia, mas **cos θ ≈ +0,04…+0,09, e −0,146 à
  superfície** — anti-alinhado, produzindo vorticidade anticiclónica.
- **Tentativa I**: deixando a tempestade ocluir livremente, o alinhamento subiu
  cerca de dez vezes e a fração streamwise foi de 0,40 para 0,64 — e só então
  |ζ| cresceu.

Ou seja, a pergunta útil não é «há cisalhamento em sentidos opostos?», mas
**«que fração dessa vorticidade está streamwise, e o gradiente do updraft está
orientado para a levantar?»**. A primeira é uma condição ambiental com muitos
falsos alarmes; a segunda foi, neste modelo, a que separou os casos.

## O instrumento

`src/storm_dynamics/tornado_precursor.py`, com duas metades deliberadamente
**não** fundidas num único número:

| função | entrada | mede |
|---|---|---|
| `environment_precursor(base)` | sondagem 1-D | ω_h, fração streamwise, SRH 0–1 e 0–3 km, cisalhamento, CAPE/CIN/LCL |
| `state_precursor(uc,vc,wc,grid)` | estado 3-D | o mesmo **mais** \|∇_h w\| e `cos θ`, sem e com condicionamento ao updraft |

A separação é o ponto: uma sondagem **não contém updraft**, logo o alinhamento
não pode ser avaliado nela, por princípio e não por omissão. O
`environment_precursor` pode dizer que um ambiente *não* produz tornado; não
pode dizer que produz.

O alinhamento é reportado também **condicionado ao updraft** (`w ≥ w_updraft`,
ponderado por w), porque a média sobre a camada dilui exatamente o ar onde a
inclinação acontece.

Uso:

```
python scripts/tornado_precursor_report.py --sounding
python scripts/tornado_precursor_report.py --sequence outputs/.../sequence.h5 --time 2790
```

## Demonstração: o ambiente diz sim, a geometria diz não

Sondagem analítica deste projeto (hodógrafo em quarto de círculo, U_max 30):

```
gates passed: omega_h_0_1km, streamwise_fraction, srh_0_1km, lcl
gates failed: none
  omega_h_0_1km_s            +0.015701
  streamwise_fraction_0_1km  +0.81153
  srh_0_1km_m2_s2            +242.42
  srh_0_3km_m2_s2            +648.19
  CAPE_J_kg                  +2225.9      LCL_m  +1068.4
```

Estado do modelo em t = 2790 s, camada 0–1000 m — o instante do pico do vórtice
de baixo nível rastreado:

```
gates passed: omega_h, streamwise_fraction
gates failed: alignment
  omega_h_mean_s        +0.015952     streamwise_fraction  +0.72546
  grad_h_w_mean_s       +7.9813e-05   alignment_mean       -0.081894
  updraft_fraction      +0.014243     alignment_updraft    -0.01586
```

O ambiente passa **todos** os portões. O estado falha no alinhamento, que é
**negativo**: naquela camada a inclinação produz vorticidade anticiclónica. Uma
previsão baseada só na regra verbal teria dito «tornado».

## Validação: o instrumento reencontra a tentativa I

Aplicado à sequência de 600 m de 10 em 10 saídas, o alinhamento condicionado ao
updraft na camada 0–1 km segue o ciclo de vida da tempestade:

| t (min) | 5 | 15 | 25 | 35 | 45 | 50 | 55 | 60 | 65 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `alignment_updraft` | −0,126 | +0,006 | −0,061 | −0,089 | −0,046 | **+0,073** | **+0,143** | +0,129 | **+0,165** |

Negativo ou nulo durante os primeiros ~45 minutos, e claramente positivo a
partir dos 50 — a fase de oclusão. Isto é a assinatura que a tentativa I
documentou por outra via («o alinhamento subiu cerca de dez vezes, 0,02 para um
pico de 0,20»), recuperada aqui por um código independente. Entretanto |ω_h|
mantém-se praticamente constante (0,0149 a 0,0160) durante todo o período: a
oferta de vorticidade horizontal nunca foi o que mudou.

Em t = 0 não há updraft e a métrica condicional devolve `nan`, como deve.

## Limites honestos

- Os limiares em `GATES` estão calibrados nas corridas **deste modelo** e
  separam os seus próprios casos. **Não são uma afirmação de skill
  climatológico**, e a amostra é pequena e de modelo único.
- `state_precursor` mede geometria num instante. Não é um prognóstico: não
  contém persistência, nem a ligação à superfície, nem intensidade. O
  classificador de seis níveis (`classification.py`) continua a ser o que
  decide se existe vórtice.
- A fração streamwise depende do quadro: passe `storm_motion` ou meça-a no
  quadro relativo ao solo, que não é o quadro que a física usa. As sequências
  deste projeto já estão no quadro relativo à tempestade, pelo que aí o valor
  correto é `(0, 0)`.
- O `environment_precursor` usa a forma ambiente `ω_h = (−∂v/∂z, ∂u/∂z)`, que
  despreza os termos em ∂w — inexistentes numa sondagem. É um limite superior do
  que a tempestade pode ingerir, não do que ela inclina.

## Artefactos

- `src/storm_dynamics/tornado_precursor.py`
- `scripts/tornado_precursor_report.py`
- `tests/test_tornado_precursor.py` — inclui a identidade de rotação
  streamwise²+crosswise² = ξ²+η², o caso de cisalhamento unidirecional (streamwise
  exatamente zero) e o caso anti-alinhado.
