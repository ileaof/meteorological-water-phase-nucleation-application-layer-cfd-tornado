# Por que o tilting muda de sinal — diagnóstico da sequência existente

Data: 6 de setembro de 2026. Escopo: pai idealizado de 600 m; principalmente 2400–3300 s, abaixo de 2 km. **Nenhuma nova simulação e nenhuma alteração do solver, da física, dos parâmetros, da LES ou da termodinâmica nesta investigação.**

**A inversão nas parcelas acompanhadas ocorre principalmente pela mudança de orientação de ∇h w, especialmente pela troca de sinal de ∂w/∂x, enquanto ξ permanece negativa.** Na parcela central, entre 2790 e 2940 s, o gradiente gira +106,06°, mas a vorticidade horizontal gira apenas +5,75°. O ângulo entre eles passa de 48,31° para 148,63°: a inclinação passa a retirar vorticidade vertical positiva. Os incrementos salvos sustentam uma combinação de gradientes de flutuabilidade, transporte da parcela através da estrutura espacial de w e ajuste de pressão. Eles não demonstram um agente físico ou numérico único responsável pela falha da tornadogênese.

O mecanismo geométrico está determinado; sua atribuição causal exclusiva e o vínculo necessário com a ausência de um tornado resolvido não estão. A resolução horizontal é 600 m, e o primeiro centro vertical está a 39,89 m: “baixo nível” não significa superfície física resolvida.

## 1. Comprovado pelos dados

### As mesmas parcelas mudam de relação geométrica com a corrente ascendente

Foram reutilizados exatamente os IDs 0–26 e as posições do diagnóstico anterior. A parcela 13 coincide, em 2790,253 s, com o máximo de ζ abaixo de 500 m do componente anteriormente rastreado. Não houve nova semeadura nem escolha de máximos independentes.

Nas 27 parcelas, o estimador principal tem a primeira amostra negativa sustentada entre **2790,253 e 2850,014 s**. Para a parcela 13, a mudança está no intervalo **2790,253–2820,426 s**. Exigiram-se duas amostras negativas consecutivas depois da última positiva entre 2700 e 2800 s. A cadência de aproximadamente 30 s não permite localizar um instante exato.

| Parcela 13 | Antes: 2700,327 s | Pico Euleriano: 2790,253 s | Após: 2940,146 s |
|---|---:|---:|---:|
| z (m) | 262,16 | 491,72 | 1066,29 |
| ξ (s⁻¹) | −0,024953 | −0,018284 | −0,010685 |
| η (s⁻¹) | +0,000222 | −0,000726 | −0,001505 |
| ∂w/∂x (s⁻¹) | −0,001369 | −0,000588 | +0,005254 |
| ∂w/∂y (s⁻¹) | −0,000389 | −0,000716 | −0,002269 |
| ξ ∂w/∂x (10⁻⁶ s⁻²) | +34,154 | +10,757 | −56,135 |
| η ∂w/∂y (10⁻⁶ s⁻²) | −0,086 | +0,520 | +3,415 |
| Tz (10⁻⁶ s⁻²) | +34,068 | +11,277 | −52,720 |
| Ângulo ωh–∇h w | 16,37° | 48,31° | 148,63° |
| w (m s⁻¹) | +1,261 | +4,122 | +1,629 |

Em 2820,426 s, ∂w/∂x já é +0,000588 s⁻¹, ξ continua −0,016595 s⁻¹ e o ângulo é 113,08°. A parcela ainda sobe a **4,579 m/s**, a 623,86 m. A inversão antecede sua entrada em movimento descendente. A parcela atinge seu próprio máximo de ζ em aproximadamente **2850 s**, depois do máximo Euleriano de 2790 s; os dois eventos não são intercambiáveis.

O termo x domina a troca de sinal. A parcela y é pequena e positiva até 2940 s, opondo-se parcialmente ao tilting negativo. η também evolui, mas sua reorganização não inicia a inversão principal.

### A decomposição geométrica identifica qual vetor domina

Usou-se T=|ωh||∇h w| cos(β−α). A decomposição simétrica de diferenças finitas separa orientação de cada vetor e módulo de cada vetor, fechando ΔT até arredondamento. É uma identidade geométrica, não um experimento de desligamento de forças.

| Contribuição a ΔT da parcela 13, 2790–2940 s | 10⁻⁶ s⁻² |
|---|---:|
| Mudança da orientação de ωh | +2,346 |
| Mudança da orientação de ∇h w | **−62,118** |
| Mudança do módulo de ωh | +2,355 |
| Mudança do módulo de ∇h w | −6,580 |
| ΔT observado | **−63,997** |

Fixando geometricamente a orientação inicial de ωh e usando a orientação final do gradiente, o cosseno seria −0,902. Fazendo somente a mudança da orientação de ωh, seria +0,737. Portanto, nesta parcela, a mudança da orientação do gradiente basta geometricamente para inverter o produto; a mudança de orientação de ωh, isoladamente, não basta.

A contribuição de orientação do gradiente supera em módulo a de ωh em **27/27 parcelas**, nas três janelas testadas: 2700–3000, 2790–2940 e 2790–3000 s. Partindo de 2700 s, a rotação isolada do gradiente produz cosseno negativo em 27/27; a rotação isolada de ωh preserva cosseno positivo em 27/27. Partindo de 2790 s, algumas parcelas já inverteram o sinal, por isso esse último contador cai para 18/27. Os 27 IDs são vizinhos correlacionados, não 27 experimentos independentes.

### Associação espacial e temporal ao ar frio

O cold pool foi diagnosticado por θv′ < −1 K no nível de **121,66 m**, comparando referências ao estado base e à média horizontal. A área com referência ao estado base cresce de **113,04 km² em 2400 s** para **169,92 km² em 2790 s**, **198,36 km² em 2940 s** e **290,52 km² em 3300 s**. Com referência à média horizontal, essas áreas são 111,60, 162,36, 189,00 e 272,88 km². Sua evolução não depende apenas de uma das referências.

A parcela 13 já apresenta θv′ < −1 K em aproximadamente **2640 s**, antes do pico. Sua projeção passa para dentro da máscara discreta do cold pool por volta de 2670 s. A diferença de um quadro resulta da máscara/transformada de distância da malha, não de uma fronteira física determinada com precisão submalha.

Na primeira amostra negativa das 27 parcelas:

- z está entre **462 e 811 m** e w entre **+3,32 e +5,59 m/s**;
- B está entre **−0,07466 e −0,05767 m/s²**;
- a contribuição térmica está entre **−0,06861 e −0,05447 m/s²**;
- a anomalia de θv no nível superficial, projetada sob cada parcela, está entre **−2,75 e −1,94 K**;
- a distância horizontal assinada à máscara é aproximadamente **−1,62 a −0,99 km**.

Assim, a inversão acontece em ar negativamente flutuante, sobre o cold pool, **ainda ascendente**. Estar sobre a máscara superficial não significa estar dentro de uma fronteira material tridimensional do cold pool. A passagem para downdraft ocorre depois: a parcela central tem w≈−1,22 m/s em 3000 s e −3,51 m/s em 3090 s.

Na parcela 13 em 2820 s, B≈−0,06637 m/s², sendo −0,06187 m/s² térmico. A carga de condensado contribui, mas não explica sozinha a flutuabilidade negativa. Os gradientes locais são Bx≈+2,76×10⁻⁵ e By≈−2,49×10⁻⁵ s⁻². Em 2940 s eles passam a aproximadamente −2,33×10⁻⁶ e −6,70×10⁻⁶ s⁻²: a interface também evolui ao longo da trajetória.

### O componente e a parcela não têm a mesma série

O componente anterior foi preservado. Seu máximo baixo de ζ muda de célula ao longo do tempo. Nesse máximo, T é negativo em 2850–2880 s, mas volta a positivo em 2910–3030 s, enquanto ζ máxima continua enfraquecendo. A parcela 13 mantém T negativo nesse intervalo. Portanto, não está comprovada uma inversão uniforme e permanente de todo o vórtice após 2790 s.

Também não existe uma única transição em toda a história: a parcela 13 tem T negativo perto de 2400–2490 s, positivo durante a aproximação ao pico, negativo na fase principal de enfraquecimento e volta a positivo por volta de 3120 s, ficando fracamente positivo em 3300 s. A pergunta sobre “após o pico” precisa manter esse recorte temporal.

## 2. Fortemente suportado

### Reorganização do gradiente por flutuabilidade e transporte, com ajuste de pressão

Foi calculado um orçamento independente de g=∇h w ao longo das posições já salvas:

`Δg = Σ interpolação[∇h(Δw por operador)] + ∫(u·∇)g dt + residual de amostragem`.

Os Δw são os incrementos nativos acumulados **em todos os passos** entre duas saídas. O transporte material foi separado em horizontal e vertical. As contribuições angulares usam `(g_x Δg_y − g_y Δg_x)/|g_médio|²`, com residual explícito da linearização em relação à rotação observada.

| Rotação de ∇h w, parcela 13 | 2790–2820 s | 2790–2940 s |
|---|---:|---:|
| Flutuabilidade | +50,30° | +55,17° |
| Advecção do solver | −20,39° | +1,03° |
| Transporte horizontal da parcela | +36,99° | +50,49° |
| Transporte vertical da parcela | −1,71° | −11,78° |
| Projeção de pressão | +4,16° | +7,54° |
| LES | −3,21° | −3,57° |
| Arrasto superficial direto | 0° | 0° |
| Residual de amostragem | +7,43° | +10,67° |
| Residual da linearização angular | −6,38° | −3,49° |
| Rotação observada | **+67,19°** | **+106,06°** |

Advecção Euleriana e transporte da parcela não são forças independentes. Sua soma representa a contribuição líquida associada a advecção/deformação nessa avaliação: **+14,89° no cruzamento** e **+39,75° até 2940 s**. Interpretar só uma parcela grande ignoraria cancelamentos. A pressão é a correção do operador de projeção; não foi separada em pressão dinamicamente forçada por rotação/deformação e pressão induzida por flutuabilidade.

Em coordenadas cartesianas, de 2790 a 2940 s, Δwx observado é aproximadamente **+0,005842 s⁻¹**. As contribuições são +0,002793 da flutuabilidade, +0,001477 da projeção, +0,001859 do transporte horizontal, +0,000366 do vertical, −0,000820 da advecção do solver, −0,000098 da LES e +0,000265 de residual. Esse fechamento com residual mostra que a conclusão não depende apenas de um ângulo normalizado.

No conjunto de 27 parcelas, a rotação angular atribuída à flutuabilidade nessa janela tem mediana +55,17°, mas varia de −12,09° a +113,83°. O transporte horizontal varia de +9,02° a +61,86°. A dominância geométrica do gradiente é uniforme; a repartição entre os agentes que o reorganizam **não é uniforme**.

Isso sustenta que as parcelas atravessam um campo ascendente em evolução, no qual a aceleração vertical diferencial e a geometria espacial mudam o sinal de wx. O episódio não exige uma reversão prévia do vetor de vorticidade horizontal, nem a entrada prévia em downdraft.

### Baroclinicidade horizontal: presente, mas não é uma simples reversão de ωh

No modelo reduzido, o rotacional da força vertical de flutuabilidade gera `(By, −Bx, 0)`, respeitando a interpolação centros–faces–centros usada na força nativa. Isso é a geração horizontal associada à flutuabilidade representada por estas equações; não é uma reconstrução completa do termo compressível `∇ρ × ∇p / ρ²`.

O orçamento horizontal anterior foi reutilizado. Entre 2790 e 2940 s, na parcela 13, a contribuição angular da flutuabilidade para ωh é **+9,56°**, a do estiramento vetorial +5,52°, a da LES +2,73° e a diferença do operador advectivo −10,43°. Coriolis contribui −1,04°, projeção −1,70°, dilatação +0,24°, residual +0,85° e linearização +0,01°, somando os **+5,75° observados**.

A interface de flutuabilidade participa quantitativamente dos dois vetores. O efeito angular acumulado sobre ∇h w é muito maior nesta janela. Não se deve interpretar as contribuições quase canceladas em ωh como ausência de geração baroclínica.

### LES e arrasto: limites quantitativos da atribuição

A LES contribui **−3,57°** para a rotação de ∇h w e **+2,73°** para ωh na parcela central entre 2790 e 2940 s. Nas 27 parcelas, sua contribuição à orientação de ∇h w fica entre −3,94° e −0,86°, oposta ao sentido dominante da rotação principal. Não há evidência de que a LES direta dispare essa inversão geométrica.

O arrasto direto amostrado é **zero nas 27 parcelas nessa janela**. Elas estão acima da camada de aplicação de 150 m. Isso exclui um impulso direto local de arrasto nessa etapa do orçamento, mas não seus efeitos herdados na vorticidade, na circulação ou na pressão. O diagnóstico anterior registra contribuições superficiais em fases mais precoces.

O orçamento anterior de ζ mostra LES negativa durante parte do enfraquecimento. Contribuir para a perda de ζ e provocar a troca de sinal de T são afirmações distintas; a segunda não é sustentada aqui para a LES.

## 3. Sugestivo

A sequência é compatível com uma reorganização da corrente ascendente associada ao crescimento do cold pool: as parcelas frias ainda ascendem, o lado do gradiente horizontal de w amostrado muda, T torna-se negativo, ζ das parcelas enfraquece e depois elas descem. Os mapas horizontais e cortes verticais mostram simultaneamente B, w, os vetores e as posições.

Essa cadeia fornece uma interpretação física coerente para a perda de manutenção da rotação baixa neste episódio. A qualificação “cold pool provoca a falha” seria mais forte que a evidência: a parcela já está fria antes da inversão, há contribuições advectivas e de pressão, e não foi executado um experimento causal. A evolução da pressão também pode comunicar influências de outras regiões.

## 4. Não determinado

- Não foi isolada uma causa única da ausência de um tornado de superfície. T negativo nas parcelas é um mecanismo de redução de ζ, mas não demonstra sozinho que toda manutenção da concentração baixa se torna impossível.
- Não foram separados evaporação, mistura e transporte como causas da anomalia térmica. A flutuabilidade foi decomposta em contribuição térmica, vapor e carga, sem alterar termodinâmica.
- Não foi determinado quanto da evolução do campo e das trajetórias resulta de efeitos indiretos anteriores de superfície ou LES.
- Não foi quantificada uma “difusão numérica pura”. A diferença entre o rotacional da advecção implementada e a forma material contínua contém efeitos da forma conservativa `−div(u u_i)` em um escoamento anelástico com `div(u)≠0`, além de discretização e amostragem. Não pode ser renomeada como difusão.
- Não foi estabelecida convergência espacial. Os estimadores alternativos usam a mesma malha de 600 m. A concordância de sinal não prova que um vórtice de escala tornádica está resolvido.
- Não foi vinculada causalmente a pressão dinâmica ao funil condensado. O diagnóstico anterior mostrou pressão termodinâmica de perturbação nula, `p_dyn` separado e sem acoplamento à saturação. O déficit local no pico era cerca de −4,70 Pa, com base de nuvem do componente em aproximadamente 706 m. Estes fatos permanecem distintos da geometria de T.

## Método, validação e arquivos

Entrada: [sequência documentada](DIAGNOSTIC_SEQUENCE.md) e [diagnóstico anterior](DIAGNOSTIC_SEQUENCE_FINDINGS.md), arquivo `outputs/diagnostic_sequence_20260905/sequence.h5`, 131 saídas de 0–3900 s. Nesta análise foram lidos 33 quadros de 2370–3330 s, incluindo margens temporais. Usou-se a camada de 0–2 km com halo já arquivada. As posições são as mesmas trajetórias RK4 anteriores, com velocidade nativa nas faces, interpolação trilinear espacial e linear temporal.

As derivadas usam diferenças centradas de segunda ordem e coordenadas z não uniformes, calculadas independentemente do operador do solver. Compararam-se três estimadores de T: produto dos vetores interpolados, interpolação do produto na malha e derivada analítica do interpolante trilinear das velocidades centradas. A concordância de sinal com o primeiro é **99,04%** no segundo e **93,55%** no terceiro, em 2400–3300 s. O terceiro desloca as primeiras amostras negativas do conjunto para 2730–2880 s; na parcela central os três preservam o bracket 2790–2820 s. RMS da diferença de T: 1,77×10⁻⁶ e 1,17×10⁻⁵ s⁻², respectivamente. A precisão temporal individual é limitada.

O orçamento discreto Euleriano anterior fecha até arredondamento, mas o orçamento material de ζ anterior tem residual RMS de **18,1%** da mudança de ζ. Os novos orçamentos angulares não substituem esse limite. Depois de aproximadamente 3120 s, o gradiente pode ficar pequeno e girar rapidamente; a linearização angular da fase pós-pico completa fica mal condicionada, com residual muito grande. **A atribuição de fontes principal foi limitada à janela 2790–2940 s**, e os incrementos cartesianos e todos os resíduos estão salvos. O ângulo e a decomposição geométrica entre estados não dependem dessa linearização de fontes.

Quatro testes analíticos passaram: fechamento da decomposição para vetores arbitrários, rotação isolada de cada vetor e derivadas trilineares em coordenadas não uniformes. Nenhum desses testes integra uma simulação atmosférica. Metadados registram hashes do código, diagnósticos de entrada e comparação dos fontes atuais com os arquivados na sequência.

Reprodução, a partir da raiz do repositório:

```powershell
python scripts/analyze_tilting_reversal.py outputs/diagnostic_sequence_20260905/sequence.h5 --prior outputs/diagnostic_sequence_20260905/analysis --out outputs/tilting_geometry_20260906
python -m pytest tests/test_tilting_geometry.py -q
```

Resultados: [tabelas automáticas](../outputs/tilting_geometry_20260906/RELATORIO_DADOS.md), [metadados](../outputs/tilting_geometry_20260906/metadata.json), [resumo numérico](../outputs/tilting_geometry_20260906/summary.json). Os CSV preservam IDs, tempos, componentes, unidades nos nomes e resíduos.

- [Séries da parcela 13](../outputs/tilting_geometry_20260906/central_geometry.png).
- [Mesmos 27 IDs](../outputs/tilting_geometry_20260906/ensemble_geometry.png).
- [Orientações e decomposição exata](../outputs/tilting_geometry_20260906/orientation_attribution.png).
- [Mapas horizontais dos vetores e de B](../outputs/tilting_geometry_20260906/horizontal_geometry.png).
- [Cortes verticais até 2 km](../outputs/tilting_geometry_20260906/vertical_sections.png).
- [Trajetórias 3D](../outputs/tilting_geometry_20260906/trajectories_3d.png): eixos em km, escalas gráficas diferentes; setas unitárias apenas direcionais, ωh preto e ∇h w verde.
- [Contribuições angulares por operador](../outputs/tilting_geometry_20260906/angular_sources.png).
- [Evolução do cold pool](../outputs/tilting_geometry_20260906/cold_pool_evolution.png).
- [Parcela versus máximo do mesmo componente](../outputs/tilting_geometry_20260906/parcel_vs_component.png).

O trabalho termina nesta atribuição diagnóstica, sem propostas de mudanças físicas.
