# O fecho subgrid não converge junto ao solo: 600 m contra 300 m por altura

Data: 2026-09-13. Nenhuma simulação foi executada. Reanálise do par de
proveniência v4 já existente; 15 s de CPU.

## Resposta principal

Refinar dx de 600 para 300 m **remove** a contribuição do fecho LES em altura,
como se espera de um termo subgrid, e **amplifica-a** junto ao solo. A quota da
LES na vorticidade vertical cai para 0,65 do valor original acima de 1,5 km e
sobe para 1,74 no nível mais baixo. Em valor absoluto, a ζ atribuída à LES a
39,9 m quase **triplica** (×2,72) enquanto o total no mesmo nível cresce apenas
×1,40.

A camada onde o vórtice tenta ligar-se ao solo é, portanto, a única onde o
fecho **não** está a convergir. Refinar horizontalmente não a vai resolver.

Isto responde à pergunta que o teste de controle deixou em aberto
([CONTROL_PARCEL_TEST.md](CONTROL_PARCEL_TEST.md)): o contraste medido abaixo de
700 m — incremento do fecho positivo, inclinação resolvida negativa — não é um
artefacto que o refinamento horizontal dissolva.

## O par usado

Duas corridas de proveniência v4 que partem **do mesmo estado** em 2790 s e
cobrem a mesma janela até 3300 s, com **grade vertical idêntica** (20 níveis,
zc de 39,9 a 2537,1 m, primeira célula de 80 m). Só dx muda.

| | 600 m | 300 m |
|---|---|---|
| ficheiro | `outputs/vorticity_provenance_long_v4_2790_3300/` | `outputs/resolution_300m_provenance_20260909/` |
| malha horizontal | 120 × 120 | 240 × 240 |
| instantes | 18 | 18 |
| passos na janela | 1015 (dt médio 0,502 s) | 1283 (dt médio 0,397 s) |
| tempo de parede | 1978 s (33 min) | 6194 s (103 min) |

A partição v4 é aditiva e exata: a soma das nove fontes reproduz `omega_total`
com erro relativo máximo de **1,50e−14** nas 36 leituras. Isso valida a leitura,
não a física.

Máscara: cilindro de raio físico de 4,2 km centrado no vórtice rastreado de cada
corrida, nível a nível. Em metros, nunca em células.

## O que os números dizem

Quota absoluta da LES, |ζ_LES| dividido pela soma de |ζ| das nove fontes, no fim
da janela:

| z (m) | 600 m | 300 m | razão | \|ζ_LES\| 600 | \|ζ_LES\| 300 | razão | razão do total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 39,9 | 0,0301 | 0,0468 | **1,56** | 4,09e3 | 1,11e4 | **2,72** | 1,40 |
| 121,7 | 0,0550 | 0,0957 | **1,74** | 9,12e3 | 2,54e4 | **2,78** | 1,54 |
| 207,5 | 0,0683 | 0,0931 | 1,36 | 1,09e4 | 2,66e4 | 2,43 | 1,54 |
| 297,7 | 0,0628 | 0,0820 | 1,31 | 1,09e4 | 2,37e4 | 2,18 | 1,51 |
| 491,7 | 0,0713 | 0,0724 | 1,02 | 1,45e4 | 1,99e4 | 1,37 | 1,34 |
| 705,7 | 0,0692 | 0,0671 | 0,97 | 1,51e4 | 2,08e4 | 1,37 | 1,24 |
| 1068,4 | 0,0569 | 0,0443 | 0,78 | 1,42e4 | 1,88e4 | 1,32 | 1,23 |
| 1488,3 | 0,0461 | 0,0315 | 0,68 | 1,55e4 | 1,63e4 | 1,05 | 1,23 |
| 1974,4 | 0,0345 | 0,0223 | **0,65** | 1,77e4 | 1,58e4 | **0,89** | 1,23 |

A média sobre os 18 instantes dá o mesmo padrão e o mesmo ponto de inversão:
1,73 a 39,9 m, 1,38 a 121,7 m, cruzando 1,0 entre 300 e 500 m, 0,65 a 2 km.

![Quota do fecho por altura](media/storm/closure_share_by_height.png)

O efeito é específico da LES. No nível mais baixo, no mesmo refinamento, a quota
do arrasto de superfície cai (0,93), a de Coriolis cai (0,59) e a da projeção
fica igual (0,95). Só a da LES sobe.

## Porque não é o passo de tempo

O par varia dx e dt em conjunto — o relatório de 300 m já assinalava isso — e a
corrida fina dá 1,26× mais passos. Se o número de passos estivesse a inflacionar
a atribuição acumulada, empurraria no mesmo sentido em todas as alturas. Ele
inverte de sinal por volta dos 500 m, sob exatamente a mesma diferença de dt.
O dt não explica o resultado.

## Interpretação, e o que nela é hipótese

Medido: a ζ atribuída ao fecho junto ao solo cresce mais depressa do que o total
quando dx é dividido por dois.

Hipótese: a grade **vertical não foi refinada**. A primeira célula tem 80 m nas
duas corridas. Ao refinar só na horizontal, o vórtice de baixo nível concentra-se
— o pico de ζ a 121,7 m é 3,10× maior a 300 m e a largura a meia altura cai para
metade — de modo que a deformação resolvida nessa camada cresce mais depressa do
que Δ² diminui, e um fecho do tipo Smagorinsky responde contribuindo **mais**.
A camada limite continua tão mal resolvida na vertical como antes, e o fecho
paga a diferença.

Isto é consistente com o que a série de tentativas A–K já tinha medido por outra
via: a resolução vertical junto à superfície dominava a ligação ao solo por um
fator de 6,6, enquanto a magnitude do arrasto era quase irrelevante.

## O próximo teste, com previsão pré-registável

A hipótese acima é falsificável e o repositório já tem o botão: `z_faces_m` em
`GridConfig`, com validação, comprometido em 7c65f7c.

**Previsão:** refinar a vertical junto à superfície, mantendo dx = 300 m, deve
fazer a quota da LES em 39,9–207,5 m **descer**, ao contrário do que fez o
refinamento horizontal. Se subir ou ficar igual, a hipótese está errada e o
problema é da forma do fecho, não da grade.

**Configuração proposta:** mesmo estado inicial de 2790 s, mesma janela de 510 s,
dx = 300 m, faces verticais explícitas com a primeira célula em cerca de 10 m no
lugar de 80 m — aproximadamente dez níveis adicionais abaixo de 500 m — e a
mesma instrumentação de proveniência v4.

**Custo estimado**, extrapolado dos tempos reais acima e não de uma suposição:
cerca de 21% mais células e um dt reduzido pelo CFL vertical, dando da ordem de
**2,5 a 3,5 horas de GPU** contra as 1,7 h da corrida de 300 m. O ficheiro de
proveniência deve passar de 6,3 GB para cerca de 9–10 GB se a coluna de saída
acompanhar os níveis novos.

**Aviso de disco:** restam 95 GB livres em C: de 932 GB (90% ocupado). Cabe, mas
a margem é estreita e deve ser reconferida antes de autorizar.

Nada disto foi executado. Fica como proposta, conforme o escopo em vigor.

## Limitações

- Um único par determinístico; não é um ensemble.
- A janela 2790–3300 s não cobre a fase preparatória de 2370–2610 s. A pergunta
  aqui respondida é a do comportamento permanente junto ao solo, não a da
  origem da semente.
- O índice de condicionamento da proveniência cresce mais depressa na corrida
  fina (κ_T 3,65 contra 2,69 no fim da janela). A atribuição é exata como
  partição, mas a sua interpretação causal degrada-se ao longo da janela nas
  duas corridas.
- A máscara cilíndrica de 4,2 km inclui estrutura fora do núcleo; a conclusão é
  sobre a camada, não sobre o núcleo isolado.

## Artefactos

- `scripts/source_attribution_by_height.py` — a reanálise e a figura.
- `outputs/source_attribution_by_height_20260913/` — CSV completo por fonte,
  nível e instante, mais o fecho da partição.
- `docs/media/storm/closure_share_by_height.csv` — tabela condensada por nível.
