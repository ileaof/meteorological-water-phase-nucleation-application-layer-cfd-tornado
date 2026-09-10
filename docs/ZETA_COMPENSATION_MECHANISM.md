# Mecanismo da diferença de vorticidade vertical entre STRONG-EVAP e WEAK-EVAP

## Resposta direta

**STRONG-EVAP não mantém `ζ_max` maior porque um stretching mais forte pague uma “dívida” de tilting mais negativo.** O resultado inesperado vem principalmente de três fatos que a comparação anterior colocava na mesma frase, mas que são quantidades diferentes:

1. **A vantagem de `ζ` já existia antes do tilting negativo.** Em 2790,253 s, STRONG já excedia WEAK em `1,02785×10⁻⁴ s⁻¹` na parcela central e no máximo Euleriano. Entre 2790 e 3000 s essa vantagem material diminuiu para `2,65662×10⁻⁵ s⁻¹`; ela não foi criada depois da inversão.
2. **O mínimo instantâneo de `Tz` da parcela central não representa o tilting integrado do componente.** Até 2850 s, o contraste integrado STRONG−WEAK de tilting na parcela central ainda era positivo, `+1,14063×10⁻⁵ s⁻¹`, isto é, menos destrutivo em STRONG. Entre 2790 e 3000 s ele se torna negativo, `−3,97775×10⁻⁵ s⁻¹`, mas, na média Euleriana do componente rastreado, o contraste acumulado de tilting é positivo, `+2,49407×10⁻⁵ s⁻¹`.
3. **A diferença tardia está concentrada no máximo, não na circulação total.** Em 3000 s, `ζ_max` é 4,35% maior em STRONG, enquanto a circulação no disco fixo de 3 km é 0,14% menor, a circulação na área do componente é somente 0,74% maior e a `ζ` média volumétrica difere apenas 0,047%. O máximo de STRONG está espacialmente mais concentrado e ocupa células onde o tilting local é relativamente mais favorável.

O orçamento material não permite atribuir a diferença remanescente em 3000 s a um termo físico dominante: o residual do contraste é `+2,65827×10⁻⁵ s⁻¹`, praticamente igual ao próprio contraste final, `+2,65662×10⁻⁵ s⁻¹`. O mecanismo resolvido e robusto é, portanto, **vantagem anterior seguida de redistribuição/concentração espacial do máximo e amostragem de uma região com tilting mais favorável**. A identificação do operador que causa toda essa redistribuição fica limitada pelo residual material e pela diferença entre a advecção discreta e sua reconstrução contínua.

## Dados e equações

Nenhuma simulação foi executada nesta etapa. Foram lidos exclusivamente os HDF5 e CSVs completos de WEAK-EVAP, CONTROL e STRONG-EVAP. Os três casos têm a mesma malha de 600 m, os mesmos 1985 passos entre 2370,212 e 3300,023 s e 32 saídas sincronizadas.

O orçamento Euleriano foi reconstruído dos incrementos nativos:

`∂ζ/∂t = transporte + stretching + tilting + dilatation + curl(LES) + curl(drag) + curl(Coriolis) + curl(projection) + curl(buoyancy) + outros + resto discreto`.

O orçamento material segue as trajetórias RK4:

`Dζ/Dt = stretching + tilting + dilatation + LES + drag + Coriolis + projection + buoyancy + resto do operador advectivo + residual material`.

O “resto do operador advectivo” é o rotacional do incremento advectivo nativo menos transporte, stretching, tilting e dilatation reconstruídos na forma contínua. Ele inclui diferenças de forma conservativa, anelasticidade, recentramento, discretização e divisão de operadores. **Ele não é sinônimo de transporte físico nem de difusão numérica.** O residual material inclui erro de trajetória, interpolação e quadratura na cadência de aproximadamente 30 s.

Os contrastes `SW = STRONG−WEAK`, `SC = STRONG−CONTROL` e `CW = CONTROL−WEAK` foram calculados para todos os termos, intervalos e 27 IDs. A tabela principal abaixo apresenta SW; os três contrastes completos estão em `lagrangian_window_contrasts.csv`, `lagrangian_phase_contrasts.csv` e `eulerian_component_contrasts.csv`.

## Orçamento material acumulado da parcela central

Valores em `10⁻⁶ s⁻¹`. “Mudança observada” é `[(ζS(t)−ζS(2790))−(ζW(t)−ζW(2790))]`; não inclui a diferença que já existia em 2790 s.

| Termo STRONG−WEAK | 2790–2850 s | 2790–2940 s | 2790–3000 s |
|---|---:|---:|---:|
| Tilting | +11,406 | −44,032 | −39,778 |
| Stretching | +21,683 | −8,074 | −27,410 |
| Dilatation | −18,128 | −27,898 | −21,008 |
| LES | −17,622 | −14,759 | −14,879 |
| Transport/advection remainder | +25,467 | +36,746 | +1,401 |
| Projection/pressure | ~0 | ~0 | ~0 |
| Buoyancy-related, direto em `ζ` | 0 | 0 | 0 |
| Drag direto | 0 | 0 | 0 |
| Coriolis | −0,315 | −0,850 | −1,128 |
| Outros operadores | 0 | 0 | 0 |
| Residual material | −15,059 | +12,774 | +26,583 |
| **Δζ observado, mudança desde 2790** | **+7,432** | **−46,092** | **−76,219** |
| Δζ já presente em 2790 | +102,785 | +102,785 | +102,785 |
| **Δζ total no fim da janela** | **+110,217** | **+56,693** | **+26,566** |

Até 2850 s, stretching, tilting e o resto advectivo favorecem STRONG, enquanto dilatation e LES se opõem. Entre 2850 e 2940 s, tilting e stretching passam a reduzir a vantagem; o resto advectivo e o residual limitam parcialmente essa redução. Até 3000 s, não existe termo físico positivo dominante: tilting, stretching, dilatation, LES e Coriolis somam contribuições desfavoráveis, e a vantagem inicial é quase totalmente consumida.

## Fases A–D

Também em `10⁻⁶ s⁻¹`, para a parcela central:

| Fase | Δζ observado | Stretching | Tilting | Dilatation | LES | resto advectivo | residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| A, 2790–2820 | +9,684 | +16,124 | +14,034 | −10,613 | −14,838 | +19,194 | −14,079 |
| B, 2820–2850 | −2,252 | +5,559 | −2,628 | −7,515 | −2,783 | +6,272 | −0,980 |
| C, 2850–2940 | −53,524 | −29,757 | −55,438 | −9,770 | +2,862 | +11,279 | +27,834 |
| D, 2940–3000 | −30,127 | −19,336 | +4,254 | +6,890 | −0,120 | −35,345 | +13,808 |

O valor mínimo posterior de `Tz` não descreve sua exposição total. STRONG tem tilting integrado menos destrutivo na fase A, quase neutro em relação a WEAK na fase B, mais destrutivo na fase C e ligeiramente mais favorável na fase D. A integral 2790–3000 ainda é mais negativa em STRONG, mas não de forma uniforme no tempo ou no espaço.

## Orçamento Euleriano do componente

Na média volumétrica do componente `ζ ≥ 0,003 s⁻¹`, os incrementos nativos fecham sem residual material. Para 2790–3000 s, o contraste acumulado SW, em `10⁻⁶ s⁻¹`, é:

| Parcela do orçamento Euleriano | SW acumulado |
|---|---:|
| mudança local observada | −9,120 |
| stretching | +4,586 |
| tilting | +24,941 |
| dilatation | −2,505 |
| transporte contínuo | −18,122 |
| resto do operador advectivo | −10,464 |
| LES | −7,521 |
| drag | +0,002 |
| Coriolis | −0,036 |
| projection/pressure | ~0 |
| buoyancy direto em `ζ` | 0 |
| fechamento | `−1,3×10⁻¹⁴` na unidade da tabela |

Essa média mostra por que o mínimo da parcela central não pode ser extrapolado para o vórtice inteiro: **no componente, o tilting é mais favorável em STRONG**, e o stretching também é levemente favorável. Transporte contínuo, resto advectivo e LES se opõem. O orçamento da célula que contém o máximo mostra o mesmo tipo de heterogeneidade, com contraste de tilting fortemente positivo depois de 2910 s. Como a célula do máximo pode mudar entre saídas, somar seus termos não constitui uma trajetória material; serve apenas para localizar onde o máximo está sendo mantido.

## Stretching, convergência e continuidade

A hipótese de compensação por `S_z = ζ ∂w/∂z` não é sustentada como mecanismo principal.

- Em 2790 s, a parcela central de STRONG tem `ζ` maior, mas `∂w/∂z` ligeiramente menor: `0,006976` contra `0,006989 s⁻¹` em WEAK.
- A partir de 2820 s, `∂w/∂z_STRONG−∂w/∂z_WEAK < 0` em 27/27 parcelas em todas as saídas, exceto seis sinais positivos dispersos já em 2790 s.
- A convergência horizontal em STRONG é menor nas 27/27 parcelas em 2790, 2820 e 2850 s; permanece menor em 26–27 das 27 até 3000 s.
- O contraste acumulado de stretching é positivo em 25/27 trajetórias até 2850 s, mas sua mediana torna-se `−2,741×10⁻⁵ s⁻¹` até 3000 s, positiva em apenas 10/27.

A decomposição exata com WEAK como referência,

`ΔS = (Δζ)(∂w/∂z)_W + ζ_W Δ(∂w/∂z) + (Δζ)Δ(∂w/∂z)`,

mostra que o primeiro termo, associado à `ζ` já maior, ajuda inicialmente, enquanto `ζ_W Δ(∂w/∂z)` é negativo e passa a dominar. Portanto, o stretching inicial é uma consequência parcial da vantagem prévia de `ζ`, não evidência de convergência mais forte causada pelo cold pool.

A divergência tridimensional de velocidade não é a restrição prognóstica relevante isoladamente no sistema anelástico. A divergência de massa `∇·(ρ₀u)` amostrada nas trajetórias permanece próxima do arredondamento, tipicamente `10⁻¹⁸–10⁻¹⁷ kg m⁻³ s⁻¹`, mesmo quando `∇·u` e a convergência horizontal não são zero.

## Transporte de `ζ`

Há uma assinatura temporária compatível com redistribuição advectiva, mas ela não permite identificar transporte como mecanismo dominante.

- Na fase C, os contrastes instantâneos de transporte horizontal são positivos em 27/27 parcelas em 2910 e 2940 s; o transporte vertical é positivo em 22/27 e 25/27, respectivamente. Isso significa menor exportação ou maior importação local em STRONG para essas amostras.
- O resto advectivo acumulado SW até 2940 s tem mediana `+3,675×10⁻⁵ s⁻¹` e é positivo em 18/27 trajetórias.
- Até 3000 s, sua mediana cai para `+7,881×10⁻⁶ s⁻¹`, positiva em 16/27, enquanto na parcela central é apenas `+1,401×10⁻⁶ s⁻¹`.
- Na média Euleriana do componente, o transporte contínuo e o resto advectivo são ambos negativos no contraste acumulado até 3000 s.

As diferenças de sinal decorrem de volumes Eulerianos, trajetórias e células do máximo amostrarem partes diferentes de um campo dipolar. Os mapas mostram lóbulos positivos e negativos adjacentes, consistentes com deslocamento/reorganização. Como o resto do operador não é transporte físico puro e o residual material é comparável ao sinal procurado, a importação adicional é **sugestiva durante 2850–2940 s**, não comprovada como explicação total.

## Geometria e circulação

| t (s) | Δ`ζ_max` SW (10⁻⁶ s⁻¹) | ΔΓ no componente (m² s⁻¹) | ΔΓ no disco de 3 km (m² s⁻¹) | Δárea a ~500 m (km²) |
|---:|---:|---:|---:|---:|
| 2790 | +102,8 | +553,6 | +185,8 | 0,00 |
| 2850 | +109,9 | +1631,8 | +1190,1 | +0,36 |
| 2940 | +136,2 | −1760,7 | −793,5 | −0,72 |
| 3000 | +254,7 | +465,5 | −62,9 | 0,00 |

Em 3000 s, o contraste da circulação em área fixa é praticamente nulo, enquanto o contraste do máximo mais que dobra em relação a 2790 s. A razão entre `ζ_max` e a `ζ` média na área do componente é 1,378 em WEAK e 1,427 em STRONG. A `ζ` média volumétrica é `0,0057557` e `0,0057584 s⁻¹`. Isso demonstra que a ordenação de `ζ_max` é principalmente uma diferença de concentração/localização, não de circulação global.

A área, o volume e a inclinação do eixo diferem pouco e nem sempre de forma monotônica. Não há evidência de um tubo inteiro sistematicamente mais intenso em STRONG. Os máximos de `w` e de convergência também permanecem vários quilômetros afastados da célula de `ζ_max`, reforçando que o máximo local não resume a geometria do updraft.

## Gradientes de flutuabilidade e vorticidade horizontal

STRONG produz inicialmente gradientes de flutuabilidade e `|ωh|` ligeiramente maiores nas trajetórias:

- em 2790 s, o contraste mediano de `|curl(B k)|` é `+9,16×10⁻⁷ s⁻²` e o de `|ωh|`, `+3,08×10⁻⁴ s⁻¹`, ambos positivos em 27/27 parcelas;
- em 2850–2940 s, o contraste de `|curl(B k)|` fica espacialmente misto;
- em 3000 s, volta a ser positivo em 27/27 parcelas.

No modelo, `curl(B k) = (∂B/∂y, −∂B/∂x, 0)`: seu termo vertical direto é exatamente zero. Essa vorticidade horizontal adicional somente influencia `ζ` depois de ser inclinada. Nas trajetórias centrais, o tilting integrado é mais destrutivo em STRONG após 2850 s, portanto não há evidência de que a geração de `ωh` pague a diferença positiva de `ζ` antes de 3000 s. Ela pode contribuir para a heterogeneidade espacial do tilting observada no componente, mas os dados não carregam rótulos de origem da vorticidade que permitam fechar essa mediação.

## Pressão e projeção

`p_dyn` e seus gradientes mudam entre os ramos; por exemplo, o gradiente amostrado na parcela central é maior em STRONG até 2850 s e volta a ser ligeiramente maior em 2940–3000 s. Entretanto, o rotacional efetivo do incremento de projeção é de ordem `10⁻¹⁸ s⁻¹` acumulada no orçamento material e fecha próximo do arredondamento no orçamento Euleriano.

Assim, a projeção não gera diretamente o contraste de `ζ`. A pressão pode reorganizar convergência e velocidade em passos posteriores, mas esse efeito indireto não pode ser separado contabilmente do estado evoluído usando apenas o curl instantâneo. `p_dyn` permaneceu exclusivamente dinâmico; não houve acoplamento à saturação.

## LES, superfície e Coriolis

- **LES:** opõe-se à diferença positiva. O contraste acumulado é negativo em 27/27 trajetórias nas três janelas principais; a mediana até 3000 s é `−1,644×10⁻⁵ s⁻¹`. LES não paga a dívida.
- **Drag:** a contribuição vertical direta é zero nas 27 trajetórias durante 2790–3000 s. Na média Euleriana do componente, o contraste acumulado é `1,5×10⁻⁹ s⁻¹`, desprezível. Não há base para atribuir um efeito indireto herdado ao drag.
- **Coriolis:** o contraste é negativo e pequeno, com mediana `−1,11×10⁻⁶ s⁻¹` até 3000 s.
- **Dilatation:** também tende a opor-se à diferença positiva, sobretudo até 2940 s.

## As 27 trajetórias

As 27 parcelas são vizinhas correlacionadas, não 27 experimentos independentes. Mesmo assim, a distribuição espacial separa sinais locais de padrões do conjunto.

| Contraste SW, mediana (10⁻⁶ s⁻¹) | 2790–2850 | 2790–2940 | 2790–3000 |
|---|---:|---:|---:|
| Δζ inicial | +92,60, 27/27 positivas | +92,60, 27/27 | +92,60, 27/27 |
| Δζ no fim | +114,53, 27/27 | +51,84, 20/27 | +0,55, 14/27 |
| mudança observada | +2,19, 15/27 | −50,54, 0/27 | −88,62, 0/27 |
| stretching | +21,33, 25/27 | −8,07, 13/27 | −27,41, 10/27 |
| tilting | +11,41, 18/27 | −44,03, 9/27 | −48,58, 9/27 |
| LES | −17,88, 0/27 | −16,37, 0/27 | −16,44, 0/27 |
| resto advectivo | +19,70, 24/27 | +36,75, 18/27 | +7,88, 16/27 |
| residual | −1,83, 10/27 | +20,50, 24/27 | +26,53, 21/27 |

Até 2850 s, a vantagem está espacialmente consistente. Em 3000 s, a mediana do contraste material é quase zero e apenas 14/27 parcelas ainda têm `ζ_STRONG > ζ_WEAK`, embora o máximo Euleriano permaneça maior. Isso é outra prova de que o resultado tardio pertence à cauda espacial/concentração, não ao conjunto material inteiro.

## Contrafactuais contábeis

O contrafactual abaixo remove um contraste por vez da soma da parcela central em 2790–3000 s. Ele não reintegra a dinâmica.

| Contraste removido | Δζ final diagnóstico sem o termo (10⁻⁶ s⁻¹) |
|---|---:|
| nenhum, observado | +26,566 |
| tilting | +66,344 |
| stretching | +53,976 |
| dilatation | +47,574 |
| LES | +41,445 |
| resto advectivo | +25,165 |
| projection/pressure | +26,566 |
| buoyancy direto | +26,566 |
| drag | +26,566 |
| Coriolis | +27,694 |
| residual | −0,017 |

Remover tilting, stretching, dilatation ou LES aumentaria contabilmente a vantagem porque seus contrastes são negativos. Remover o resto advectivo quase não muda o resultado em 3000 s. Remover o residual elimina toda a vantagem remanescente. Isso impede classificar qualquer termo físico como compensação dominante na parcela central até 3000 s.

## Controle numérico

- O fechamento máximo do orçamento Euleriano médio no componente é `8,52×10⁻²⁰ s⁻¹` por intervalo, muito menor que os contrastes físicos. Esse fechamento verifica que nenhum incremento de velocidade foi omitido.
- O orçamento material é menos preciso: em 3000 s, o contraste residual central de `+26,583×10⁻⁶ s⁻¹` é do tamanho do contraste final. Na mediana das 27 parcelas, o residual é `+26,526×10⁻⁶ s⁻¹`, enquanto o contraste final mediano é apenas `+0,545×10⁻⁶ s⁻¹`.
- Reduzir o passo RK4 de 2 para 1 s mudou as trajetórias em menos de 0,031 m; usar cadência de 60 s em vez de 30 s gerou medianas de 6,3–6,4 m e máximos abaixo de 19,7 m. A posição é estável, mas a quadratura de termos grandes que se cancelam continua sensível à cadência.
- A sequência de `dt` de CONTROL foi imposta aos três ramos. Ela excedeu o limite adaptativo próprio de WEAK em no máximo 0,0924% e de STRONG em 0,2304%; nenhum passo excedeu 1%. Essa assimetria é pequena, mas não existe um ramo contrafactual com `dt` adaptativo próprio para quantificar seu efeito sobre sinais de ordem `10⁻⁵ s⁻¹`. Portanto, não se pode excluí-la na escala do residual material.
- A malha de 600 m limita posição, área e circulação a estruturas resolvidas. A circulação no disco de 3 km usa dezenas de células e é mais estável que uma única célula, mas ainda é uma integral de `ζ` reconstruída na malha.

## Respostas às dez perguntas

1. **A compensação vem principalmente de stretching?** Não. Ele ajuda até 2850 s por causa da `ζ` inicial maior, mas passa a ser desfavorável; convergência e `∂w/∂z` são menores em STRONG.
2. **Existe maior convergência em STRONG?** Não nas trajetórias equivalentes. O contraste é negativo quase unanimemente.
3. **Existe transporte adicional de `ζ`?** Há uma assinatura local temporária em 2850–2940 s, mas ela não é persistente nem separável do resto discreto e do residual. Classificação: sugestivo.
4. **O tilting integrado realmente é mais destrutivo em STRONG?** Até 2850 s, não; até 2940 e 3000 s, sim na parcela central e na mediana das trajetórias. Na média do componente Euleriano, o contraste tem sinal oposto e é mais favorável em STRONG.
5. **LES contribui significativamente para o contraste?** Sim como oposição: remove relativamente mais `ζ` em STRONG. Não explica a vantagem positiva.
6. **A pressão/projeção contribui?** O curl direto não contribui acima do arredondamento. Efeitos indiretos da pressão não são isoláveis neste orçamento.
7. **Existe aumento de circulação ou apenas maior concentração?** Em 3000 s, a circulação em área fixa é praticamente igual e o máximo é 4,35% maior. Predomina concentração/localização do máximo.
8. **O mecanismo aparece nas 27 trajetórias?** A vantagem inicial aparece em 27/27. Em 3000 s, somente 14/27 mantêm sinal positivo; a concentração do máximo não representa todo o conjunto.
9. **O contraste é maior que as incertezas numéricas?** O contraste do máximo Euleriano é muito maior que o fechamento Euleriano. A atribuição material termo a termo não excede com segurança o residual acumulado em 3000 s.
10. **Existe evidência suficiente para identificar um mecanismo dominante?** Existe evidência suficiente para identificar vantagem prévia e concentração espacial como explicação do `ζ_max` maior. Não existe evidência suficiente para atribuir essa concentração a um único operador físico ou numérico.

## Produtos

Os resultados estão em `outputs/zeta_compensation`:

- `lagrangian_interval_budgets.csv`, `lagrangian_phase_contrasts.csv` e `lagrangian_window_contrasts.csv`;
- `ensemble_window_summary.csv`;
- `eulerian_budget_*.csv`, `eulerian_component_contrasts.csv` e `eulerian_component_window_contrasts.csv`;
- `component_geometry_*.csv` e `parcel_samples_*.csv`;
- `stretching_factorization.csv`;
- `diagnostic_counterfactuals.csv`;
- `summary.json` e 18 figuras automáticas em `figures/`.

As figuras numeradas incluem `ζ`, `Dζ/Dt`, stretching, tilting, LES, transporte horizontal/vertical, projeção/pressão, contrastes instantâneos e acumulados, waterfall em 3000 s, convergência e `∂w/∂z`, circulação, área/volume, `ωh` e `curl(B k)`, mapas, cortes verticais, as 27 parcelas e o diagrama final.

## COMPROVADO PELOS DADOS

- STRONG já tinha `ζ` maior em 2790 s, antes da fase principal de tilting negativo.
- A vantagem material diminuiu de 2790 a 3000 s; não foi criada pelo período posterior à inversão.
- STRONG não apresenta maior convergência nem maior `∂w/∂z` nas trajetórias equivalentes.
- Stretching, LES, dilatation e Coriolis não pagam a diferença positiva acumulada até 3000 s; LES se opõe em 27/27 trajetórias.
- O curl direto de buoyancy em `ζ` é zero e o de projection/pressure é desprezível.
- Em 3000 s, a diferença de `ζ_max` não corresponde a aumento comparável de circulação ou de `ζ` média: é uma concentração espacial do máximo.
- O tilting relevante ao componente é espacialmente heterogêneo; seu contraste na região do máximo e na média Euleriana não tem o mesmo sinal do mínimo da parcela central.

## FORTEMENTE SUPORTADO

- O aparente paradoxo decorre principalmente de comparar uma vantagem anterior e um máximo espacial com o mínimo instantâneo de `Tz` em uma trajetória.
- A redistribuição do campo coloca o máximo de STRONG em uma parte relativamente mais favorável do tilting, enquanto a circulação total permanece quase igual.
- Até 2850 s, stretching e resto advectivo ajudam a preservar a vantagem inicial; depois, o stretching deixa de ajudar.
- A ordenação tardia de `ζ_max` representa a cauda espacial do componente, não uma resposta uniforme das 27 parcelas.

## SUGESTIVO

- Transporte horizontal e vertical relativamente mais favorável em partes de 2850–2940 s pode participar da concentração.
- A maior vorticidade horizontal e os gradientes de `B` em partes da janela podem contribuir para a heterogeneidade espacial do tilting.
- Mudanças em `p_dyn` podem reorganizar indiretamente o escoamento, embora seu curl direto seja nulo.

## NÃO DETERMINADO

- Qual operador individual causa toda a concentração do máximo entre 2940 e 3000 s.
- Quanto da pequena diferença material é efeito da cadência de 30 s, do resto advectivo discreto ou da sequência de `dt` imposta.
- A mediação exata `curl(B k) → ωh etiquetada por origem → tilting dessa parcela de ωh`; faltam traçadores diagnósticos passivos que conservem a origem de cada componente de vorticidade. A instrumentação mínima seria acumular, ao longo das trajetórias, `ωh` gerada por cada operador e seu tilting subsequente, sem realimentar o estado.
- O efeito indireto isolado de `p_dyn` ou drag herdado. Separá-los exigiria uma decomposição de resposta/tangente ou novos ramos controlados, não apenas campos adicionais; nenhum ramo foi executado.
- Estruturas submalha, tornado ou funil visível.

Nenhum parâmetro, solver, LES, drag, hodógrafa, resolução ou termo termodinâmico foi alterado. Nenhuma nova simulação ou matriz causal foi executada.
