# Comparação de tornadogênese: 600 m versus 300 m

Data da conclusão: 2026-09-09.

## Resposta principal

O refinamento horizontal de 600 para 300 m reduziu claramente o déficit de
concentração da vorticidade de baixo nível. Ele produziu um núcleo mais estreito
e intenso, maior circulação e maior velocidade tangencial nas mesmas escalas
físicas. O resultado não demonstra tornadogênese resolvida: o RMW ainda possui
aproximadamente três células, o vórtice enfraquece durante a janela e a
comparação inclui uma diferença numérica ligada ao passo de tempo e ao damping
superior.

A classificação anterior continua válida nos dois casos:
`MIXED: CONCENTRATION DEFICIT DOMINANT, SOURCE-INVENTORY LOSS SECONDARY`.
No caso de 300 m, porém, o déficit de concentração é menor.

## Validade numérica de 300 m

A proveniência v4 percorreu os 1.283 passos nativos entre 2790,579551 e
3300,023494 s em GPU. Os 12 campos prognósticos permaneceram bit a bit idênticos
ao controle.

- Gate: `PASS`.
- Máximo fechamento relativo RMS de omega: `6,7982e-15`.
- Máximo erro absoluto de omega: `9,8619e-16`.
- Máximo fechamento relativo RMS de tilting: `8,7053e-15`.
- Máximo erro absoluto de tilting: `7,5081e-18`.
- Índice de condicionamento máximo: omega `2,0086`; tilting `3,6488`.
- Máximo resto de advecção: `8,7386e-15` entre as quatro medidas registradas.
- Arquivo: `provenance_long_v4.h5`, 6.313.607.913 bytes,
  SHA-256 `fc8c5e3c57fad248294fa1e7c21328ae57d2fac1639eb320b14fcadfd10805b0`.

O reparo do último snapshot alterou somente o atributo `step`, de 7398 para
7397. Todos os datasets conservaram os hashes. O SHA-256 da sequência passou de
`40680e4995b55ba76d86c5fb819f5345768ea6007a4335eb357389dba9c4e736`
para `2ade1b0911b453bba1fa808dc3871a8dc4e4f8bad612effa35b13c3368ec8fff`.

## Comparação em unidades físicas

Foram usados 18 tempos comuns de 2820 a 3300,023494 s, sem extrapolação, nível
primário de 121,659 m, bins radiais de 600 m, cilindro de raio 4,2 km e volume
recortado entre 0 e 2 km.

| Mediana temporal | 600 m | 300 m | razão 300/600 | ocorrência em 300 m |
|---|---:|---:|---:|---:|
| zeta máxima a 121,7 m (s^-1) | 0,00474 | 0,01467 | 3,10 | 18/18 maior |
| zeta máxima entre 0 e 2 km (s^-1) | 0,01356 | 0,01876 | 1,38 | 18/18 maior |
| circulação assinada, raio 4,2 km (m2 s^-1) | 18.874 | 24.999 | 1,32 | 18/18 maior |
| circulação absoluta, raio 4,2 km (m2 s^-1) | 43.525 | 57.682 | 1,33 | 18/18 maior |
| integral absoluta de zeta, 0–2 km (m3 s^-1) | 3,748e8 | 4,889e8 | 1,30 | 18/18 maior |
| integral assinada de zeta, 0–2 km (m3 s^-1) | 1,417e8 | 1,289e8 | 0,91 | 3/18 maior |
| velocidade tangencial máxima (m s^-1) | 1,92 | 2,52 | 1,31 | 18/18 maior |
| largura de zeta à meia altura (m) | 2.144 | 1.148 | 0,54 | 0/18 maior |
| RMW (m) | 919 | 893 | 0,97 | 3/18 maior |
| raio equivalente do núcleo (m) | 1.265 | 1.233 | 0,97 | 15/18 maior |
| largura da convergência à meia altura (m) | 5.288 | 4.504 | 0,85 | 0/18 maior |
| largura da queda de pressão à meia altura (m) | 4.617 | 4.537 | 0,98 | 4/18 maior |
| stretching condicional médio (s^-2) | 2,710e-5 | 3,447e-5 | 1,27 | 18/18 maior |

O aumento simultâneo dos picos, da circulação assinada e absoluta e da
velocidade tangencial, junto da redução física da largura de zeta, comprova a
melhor concentração. A integral assinada de zeta em todo o volume 0–2 km cai
9%, enquanto a integral absoluta sobe 30%. Portanto, o refinamento também
resolve mais estrutura de sinais opostos; ele não cria simplesmente mais
inventário líquido no cilindro.

## Adequação da resolução

As quantidades em células são usadas apenas para avaliar se a estrutura está
resolvida. Elas foram retiradas da tabela de mudança física.

- RMW mediano: 1,53 células em 600 m e 2,98 células em 300 m.
- diâmetro equivalente do núcleo: 4,22 e 8,22 células.
- largura de zeta à meia altura: 3,57 e 3,83 células.
- largura da convergência à meia altura: 8,81 e 15,01 células.

O RMW melhora, mas cerca de três células ainda é marginal. A largura física de
zeta diminui quase pela metade e permanece com aproximadamente quatro células;
isso impede declarar convergência de malha. Um próximo teste de concentração
deveria usar resolução horizontal menor que 300 m, sem mudar a física.

## Fontes, conversão e perda

No orçamento Euleriano de 300 m, as médias na máscara móvel são: stretching
`+2,879e-5`, tilting `+2,091e-5`, transporte `-3,997e-5`, dilatação
`-4,047e-6`, LES `-7,939e-7`, arrasto `-4,039e-8` e Coriolis
`+3,577e-7 s^-2`. Stretching e tilting estão disponíveis; transporte e
dilatação compensam sua produção, e LES é um sumidouro secundário.

Em 3300 s, 58,34% da norma positiva de zeta de 300 m pertence ao estado
antecedente `initial` e 37,66% à projeção. Na vorticidade horizontal absoluta,
37,93% é antecedente, 33,09% é projeção e 23,48% vem de buoyancy. A contribuição
vertical direta de buoyancy é zero pela geometria do termo; sua contribuição
aparece primeiro na vorticidade horizontal e pode alcançar zeta por tilting.
O rótulo `projection` mede a redistribuição pelo incremento de pressão e não
deve ser interpretado isoladamente como uma fonte ambiental independente.

Os dois casos perdem força na janela. A razão final/pico de zeta é 0,703 em
600 m e 0,695 em 300 m; para a circulação de 4,2 km é 0,440 e 0,370. Assim, a
perda de inventário continua como mecanismo secundário mesmo após melhorar a
concentração.

## Coerência vertical e domínio

A inclinação mediana da superfície a 2 km cai de 1,72 para 1,51 km. O máximo
deslocamento do eixo em qualquer altura cai de 4,42 para 1,59 km. Há melhora de
coerência no caso refinado, mas as trajetórias não são idênticas e o máximo de
300 m chegou a 7,29 km no começo da janela.

O domínio simulado em ambos os casos é 72 x 72 x 15 km. O cold pool de -1 K
permanece pelo menos 11,1–11,25 km das bordas, mas condensado toca a costura
periódica e alcança a camada superior. Updrafts de 5 m s^-1 chegam a 0,3 km da
borda em 600 m e 1,05 km em 300 m. Isso não prova auto-interação, mas impede
certificar que o domínio seja irrelevante para a supercélula completa.

Um domínio de 120 x 120 x 20 km adicionaria cerca de 24 km de margem em cada
lado e 5 km no topo, sendo o próximo teste apropriado para extensão do domínio.
Ele é provavelmente suficiente para esta tempestade idealizada durante a janela
analisada, mas isso permanece uma hipótese até uma comparação de domínio. A
região superior começa em 12,74 km e é alcançada pela nuvem. A condição de
damping é aplicada quatro vezes por passo, sem escala por dt; a exposição
nominal é 27% maior em 300 m. Por isso a experiência mede conjuntamente
resolução e efeitos numéricos dependentes do passo de tempo.

## Força das conclusões

- **Comprovado pelos dados:** o gate numérico passou; 300 m concentra melhor o
  vórtice em unidades físicas; o inventário absoluto e a circulação aumentam;
  ambos os casos enfraquecem; o condensado alcança bordas laterais e o topo.
- **Fortemente suportado:** a deficiência dominante em 600 m era concentração e
  foi parcialmente reduzida em 300 m; transporte/dilatação e LES contribuem para
  o enfraquecimento; a coerência vertical melhora.
- **Sugestivo:** 120 x 120 x 20 km reduzirá a influência geométrica das bordas e
  do topo; resolução abaixo de 300 m poderá concentrar ainda mais o núcleo.
- **Não determinado:** convergência de malha, irrelevância causal do damping,
  suficiência definitiva de 120 x 120 x 20 km e tornadogênese resolvida.

## Artefatos

- Proveniência: `outputs/resolution_300m_provenance_20260909/summary.json` e
  `provenance_long_v4.h5`.
- Diagnóstico de 300 m: `outputs/resolution_300m_audit_20260909/`.
- Comparação corrigida: `outputs/resolution_comparison_20260909/summary.json`,
  `metric_comparison.csv`, `resolution_adequacy_cells.csv` e `manifest.json`.
- Auditorias de domínio: `outputs/domain_extent_audit/` e
  `outputs/domain_extent_300m_audit/`.
