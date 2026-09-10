# Experimento causal piloto da evaporação da chuva

## Resposta direta

**Parcialmente.** A modificação controlada da evaporação da chuva produziu uma alteração causal, mensurável e ordenada na evaporação efetivamente realizada e nas métricas areais e médias do cold pool. Ela também produziu diferenças pequenas em `∇h w` e na magnitude de `Tz`. Contudo, **não alterou o instante resolvido da inversão sustentada de `Tz` e não reduziu a manutenção de `ζ` no sentido previsto por H1/H2**. Nos três casos, a parcela central passou de `Tz > 0` em 2790,253 s para `Tz < 0` em 2820,426 s. Depois da inversão, `STRONG-EVAP` apresentou `ζ` baixa ligeiramente maior, não menor, que `CONTROL` e `WEAK-EVAP` em todas as sete saídas entre 2820 e 3000 s.

Assim, o piloto sustenta causalmente a cadeia

`f → evaporação real → ΔT/Δqv/Δqr → Δθv/ΔB → cold pool`,

detecta uma modulação fraca da geometria de `∇h w` e de `Tz`, mas **não sustenta a cadeia causal completa até uma perda antecipada ou mais intensa de `ζ` de baixo nível**.

## Escopo do teste

Este é um experimento causal de reinício, específico da fase crítica. Os três ramos partem do mesmo estado arquivado e bit a bit idêntico em `t = 2370,211552 s`; o fator passa a atuar somente a partir desse instante e os ramos terminam em `t = 3300,023494 s`. Cada ramo contém 1985 passos e 32 estados tridimensionais sincronizados, com cadência próxima de 30 s. O teste estima a resposta marginal da tempestade já formada durante a janela da inversão. Ele não testa o efeito do fator desde a iniciação da convecção.

Os casos congelados após a triagem foram:

| Caso | `rain_evaporation_factor` |
|---|---:|
| WEAK-EVAP | 0,95 |
| CONTROL | 1,00 |
| STRONG-EVAP | 1,05 |

A triagem avaliou `δ = 0,05`, `0,10` e `0,20`. Em 2790 s, `δ = 0,05` já separava WEAK e STRONG por 11,48 milhões de kg de chuva evaporada e por 2,88 km², ou oito células horizontais, na área com `θv′ < −1 K`. A diferença máxima de `w` era 0,0075% do controle e a de condensado total, 0,0267%. Como `δ = 0,05` era o menor candidato que atendia ao critério pré-definido, os candidatos maiores foram descartados e não houve segunda escolha depois de observar `ζ`.

Depois de remover o próprio fator, as configurações serializadas dos três ramos são idênticas. Malha, estado inicial, hodógrafa, LES, arrasto, microfísica restante, dispositivo GPU e sequência de `dt` são os mesmos. O ramo `CONTROL` reproduziu bit a bit, com erro máximo zero, todas as 31 saídas históricas posteriores ao reinício até 3300 s.

## Como as quantidades foram medidas

- O componente ciclônico usa `ζ ≥ 0,003 s⁻¹`, conectividade de seis vizinhos e associação somente por sobreposição. A semente comum em 2400 s caiu na mesma célula nos três casos: `(48,9; 38,7; 1,488) km`.
- `ζ`, `ξ`, `η`, `∂w/∂x`, `∂w/∂y`, stretching e `Tz = ξ ∂w/∂x + η ∂w/∂y` foram reconstruídos independentemente das velocidades nativas.
- Foram integradas 27 trajetórias RK4 em torno do pico baixo de cada componente equivalente. A parcela 13 é a amostra central apresentada nas séries; a atribuição geométrica também foi conferida no conjunto de 27 parcelas.
- O cold pool foi medido no centro vertical mais próximo de 100 m, `z = 121,66 m`, com limiares `θv′ < −0,5`, `−1` e `−2 K`. `θv′` é relativo ao estado base. A profundidade exige conexão vertical contínua desde o primeiro nível e é truncada em 2 km.
- A “velocidade da borda” é a derivada do raio de área equivalente `sqrt(A/π)`. Ela não é uma velocidade local da frente de rajada.
- Stretching e demais tendências são integrais nativas de cada intervalo divididas pelo `dt` real. `Tz` das figuras geométricas é instantâneo.
- Comparações por evento usam brackets das saídas de 30 s. Instantes iguais nas tabelas significam ausência de deslocamento detectável nessa cadência, não igualdade sub-30 s.

## 1. A intervenção realmente alterou a evaporação?

**Sim.** A ordenação `WEAK-EVAP < CONTROL < STRONG-EVAP` ocorreu em todas as 32 saídas:

| Métrica acumulada em 3300 s | WEAK-EVAP | CONTROL | STRONG-EVAP |
|---|---:|---:|---:|
| chuva evaporada (kg) | 4,19856×10⁸ | 4,35856×10⁸ | 4,51404×10⁸ |
| resfriamento latente associado (J) | 1,05006×10¹⁵ | 1,09008×10¹⁵ | 1,12896×10¹⁵ |

O contraste STRONG menos WEAK foi `+3,15481×10⁷ kg`, ou 7,24% do valor de CONTROL. A resposta não precisa ser exatamente 10% porque o fator multiplica a taxa antes dos limitadores; saturação, chuva disponível e a evolução posterior retroagem sobre a taxa realizada.

No domínio de 0–2 km, em 2790 s, STRONG menos WEAK apresentou média de `ΔT = −0,00146 K`, `Δθ = −0,00150 K`, `Δqv = +4,62×10⁻⁷ kg kg⁻¹` e `Δqr = −3,89×10⁻⁷ kg kg⁻¹`. Os RMS correspondentes foram 0,00981 K, 0,00999 K, `4,67×10⁻⁶` e `2,74×10⁻⁶ kg kg⁻¹`. Os sinais médios são os esperados para mais evaporação: mais vapor, menos chuva e resfriamento.

## 2. Alterou o cold pool?

**Sim.** A área com `θv′ < −1 K` e a intensidade média da região fria ficaram ordenadas em todas as 32 saídas. No pico de `ζ` baixa, em 2790 s:

| Métrica em 2790 s | WEAK-EVAP | CONTROL | STRONG-EVAP |
|---|---:|---:|---:|
| área `θv′ < −1 K` (km²) | 168,12 | 169,92 | 171,00 |
| `θv′` médio na região fria (K) | −1,8660 | −1,8865 | −1,9099 |
| `B` médio na região fria (m s⁻²) | −0,06052 | −0,06119 | −0,06194 |
| máximo `|∇h B|` (s⁻²) | 6,945×10⁻⁵ | 7,009×10⁻⁵ | 7,072×10⁻⁵ |
| profundidade fria média (m) | 725,6 | 731,9 | 739,8 |

O maior valor observado da área de `θv′ < −1 K` foi 284,40, 290,52 e 293,76 km², respectivamente. Esses máximos ocorreram na última saída e, portanto, são censurados à direita: são máximos dentro da janela, não máximos completos do ciclo de vida. Os mínimos de superfície observados também foram ordenados: `θv′ = −3,598`, `−3,626`, `−3,662 K` e `B = −0,1171`, `−0,1179`, `−0,1191 m s⁻²`.

A profundidade máxima atingiu o teto diagnóstico de 2 km nos três casos e não discrimina os ramos. A profundidade média no momento do pico de `ζ`, que não satura no teto, mostra a resposta ordenada acima.

## 3. Alterou `∇h w`?

**Sim, mas pouco e sem ordenação persistente ao longo de toda a fase.** Em 2790 s, na parcela central equivalente, `|∇h w|` foi `0,912`, `0,927` e `0,940 ×10⁻³ s⁻¹`: STRONG excedeu WEAK em 3,09%. Entre 2790 e 2940 s, a orientação de `∇h w` girou 105,74°, 106,06° e 106,50°. A diferença total de rotação entre STRONG e WEAK foi somente 0,76°.

A ordem da magnitude de `|∇h w|` muda depois de aproximadamente 2880 s. Considerando as 32 saídas, a ordem crescente ocorreu em 19 e a decrescente em 12; portanto, o efeito não é uma intensificação geométrica monotônica durante toda a janela. No primeiro instante negativo de `Tz`, o ângulo `ωh–∇h w` foi 113,47°, 113,08° e 112,67°: a inversão geométrica ocorreu em todos os casos, com diferença menor que 1°.

## 4. Alterou `Tz`?

**Alterou ligeiramente a magnitude, mas não o instante resolvido da inversão.** A parcela central permaneceu positiva em 2790,253 s e negativa em 2820,426 s nos três ramos. Os três estimadores usados — produto dos campos derivados, campo `T` interpolado e derivada analítica do interpolante trilinear — preservaram a inversão do conjunto; a dispersão dos primeiros cruzamentos entre as 27 parcelas também se sobrepõe nos três casos.

O mínimo pós-pico de `Tz` foi:

| WEAK-EVAP | CONTROL | STRONG-EVAP |
|---:|---:|---:|
| −5,2818×10⁻⁵ s⁻² | −5,3167×10⁻⁵ s⁻² | −5,3493×10⁻⁵ s⁻² |

STRONG foi 1,28% mais negativo que WEAK no mínimo. Porém, no próprio primeiro instante negativo, STRONG era ligeiramente menos negativo (`−8,208×10⁻⁶ s⁻²`) que WEAK (`−8,244×10⁻⁶ s⁻²`). A ordem esperada por H1 apareceu em quatro das sete saídas entre 2820 e 3000 s, insuficiente para afirmar uma resposta monotônica da fase inteira.

A decomposição exata de 2790–2940 s manteve o resultado do diagnóstico anterior em todos os ramos: a parcela central girou `ωh` apenas 5,92°, 5,75° e 5,58°, enquanto `∇h w` girou cerca de 106°. No conjunto, a parcela de orientação de `∇h w` dominou a mudança de `Tz` em 27/27 trajetórias em cada caso. Logo, o fator modulou levemente **o mesmo mecanismo de inversão**, mas não o deslocou em uma saída de 30 s.

## 5. Alterou `ζ`?

**Sim, em pequena magnitude, mas no sentido oposto à perda antecipada proposta por H1.** O pico Euleriano abaixo de 500 m ocorreu em 2790,253 s nos três casos e foi:

- WEAK-EVAP: `0,0073232 s⁻¹`;
- CONTROL: `0,0073752 s⁻¹`;
- STRONG-EVAP: `0,0074260 s⁻¹`.

O pico material da parcela central ocorreu em 2850,014 s nos três casos: `0,0080895`, `0,0081456` e `0,0081997 s⁻¹`. STRONG excedeu WEAK por cerca de 1,4% em ambas as medidas. Nas sete saídas pós-inversão até 3000 s, a relação foi `ζ_STRONG > ζ_CONTROL > ζ_WEAK`, de modo que a ordenação requerida por H1 ocorreu em zero de sete saídas.

A primeira queda de pelo menos 20%, mantida por três saídas, começou em 3000,139 s em WEAK e em 3030,422 s em CONTROL e STRONG. Isso também não segue a previsão “STRONG perde antes, WEAK perde depois”. Os três componentes continuavam conectados ao primeiro nível em 3300 s; a duração observada da conexão baixa é `≥ 689,65 s`, e o instante real de perda permanece indeterminado.

## 6. Houve resposta monotônica?

Houve resposta monotônica robusta até o cold pool:

- evaporação acumulada: 32/32 saídas ordenadas;
- área fria: 32/32;
- intensidade fria média: 32/32.

Na segunda metade da cadeia, a resposta deixou de ser persistentemente monotônica:

- `Tz` no sentido de H1: 16/32 saídas no total e 4/7 entre 2820 e 3000 s;
- `ζ` no sentido de H1: 11/32 no total e 0/7 entre 2820 e 3000 s;
- momento da inversão: o mesmo bracket de 30 s nos três casos;
- momento do pico baixo e do pico material: a mesma saída nos três casos.

Portanto, existe monotonicidade causal no processo alvo e no cold pool, mas não existe a ordenação downstream exigida para sustentar H1/H2 integralmente.

## 7. A cadeia causal completa foi sustentada?

**Não.** O experimento comprova que aumentar `f` aumenta a evaporação realizada, resfria e umidifica o estado de 0–2 km, reduz `qr`, torna as métricas médias do cold pool mais negativas e aumenta sua área. Há também diferenças pequenas e resolvidas em `∇h w` e em `Tz`. A cadeia quebra como explicação da perda de rotação porque:

1. a inversão de `Tz` ocorre no mesmo bracket em todos os casos;
2. a alteração da rotação de `∇h w` é menor que 1° entre os extremos;
3. a ordenação de `Tz` não persiste por toda a fase negativa;
4. `ζ` responde no sentido oposto ao previsto, ficando um pouco maior em STRONG após a inversão;
5. WEAK, não STRONG, cruza primeiro o critério de perda persistente.

O resultado não prova ausência de qualquer influência da evaporação sobre `ζ`. Ele mostra que, para esta intervenção de ±5%, iniciada em 2370 s, **o fortalecimento mensurável do cold pool não é suficiente para explicar causalmente o horário da inversão nem a perda subsequente de `ζ` pela cadeia unidirecional proposta**.

## Marcos absolutos e alinhados por evento

| Evento | WEAK-EVAP | CONTROL | STRONG-EVAP |
|---|---:|---:|---:|
| A. início da intensificação | ≤2370,212 s, censurado | ≤2370,212 s, censurado | ≤2370,212 s, censurado |
| B. primeiro contato com o primeiro nível | 2610,372 s | 2610,372 s | 2610,372 s |
| C. pico de `ζ` abaixo de 500 m | 2790,253 s | 2790,253 s | 2790,253 s |
| D. pico material de `ζ`, parcela 13 | 2850,014 s | 2850,014 s | 2850,014 s |
| E. inversão sustentada de `Tz` | 2790,253–2820,426 s | 2790,253–2820,426 s | 2790,253–2820,426 s |
| F. início da perda persistente de `ζ` | 3000,139 s | 3030,422 s | 3030,422 s |
| G. última conexão baixa observada | ≥3300,023 s | ≥3300,023 s | ≥3300,023 s |
| H. maior área fria observada | 3300,023 s, censurado | 3300,023 s, censurado | 3300,023 s, censurado |

O alinhamento por C e E não revela um deslocamento oculto: como os eventos centrais caem nas mesmas saídas, as séries alinhadas mantêm a mesma separação pequena observada em tempo absoluto. Isso limita a conclusão temporal à cadência de saída.

## Tabela final

| Métrica | WEAK-EVAP | CONTROL | STRONG-EVAP |
|---|---:|---:|---:|
| `f` | 0,95 | 1,00 | 1,05 |
| evaporação acumulada (kg) | 4,19856×10⁸ | 4,35856×10⁸ | 4,51404×10⁸ |
| área fria máxima observada (km²) | 284,40 | 290,52 | 293,76 |
| `θv′` mínimo a 121,7 m (K) | −3,5977 | −3,6263 | −3,6617 |
| `B` mínimo a 121,7 m (m s⁻²) | −0,11708 | −0,11791 | −0,11906 |
| tempo do pico de `ζ` baixa (s) | 2790,253 | 2790,253 | 2790,253 |
| pico de `ζ` baixa (s⁻¹) | 0,007323 | 0,007375 | 0,007426 |
| tempo do pico material (s) | 2850,014 | 2850,014 | 2850,014 |
| pico material de `ζ` (s⁻¹) | 0,008089 | 0,008146 | 0,008200 |
| primeira amostra negativa de `Tz` (s) | 2820,426 | 2820,426 | 2820,426 |
| `Tz` mínimo pós-pico (s⁻²) | −5,2818×10⁻⁵ | −5,3167×10⁻⁵ | −5,3493×10⁻⁵ |
| rotação de `∇h w`, 2790–2940 (°) | 105,74 | 106,06 | 106,50 |
| rotação de `ωh`, 2790–2940 (°) | 5,92 | 5,75 | 5,58 |
| início da perda persistente de `ζ` (s) | 3000,139 | 3030,422 | 3030,422 |
| última conexão baixa | ≥3300,023 s | ≥3300,023 s | ≥3300,023 s |
| duração da conexão baixa | ≥689,65 s | ≥689,65 s | ≥689,65 s |

## Efeitos microfísicos secundários

Essas diferenças são respostas indiretas posteriores, não ação direta do parâmetro. Entre STRONG e WEAK, a maior diferença relativa acumulada fora da evaporação da chuva foi no melting de neve, 1,04% do CONTROL. Autoconversão diferiu 0,322%, acreção 0,295%, condensação 0,214%, melting de graupel 0,191%, evaporação de nuvem 0,116% e freezing de chuva 0,091%. A precipitação superficial acumulada de chuva foi `1,04294×10⁹`, `1,03573×10⁹` e `1,02875×10⁹ kg`, uma redução de 1,37% entre WEAK e STRONG, coerente com a alteração de `qr` e com a sedimentação posterior.

Processos de gelo não receberam o fator. As diferenças neles surgem porque `T`, `qv`, `qr` e o escoamento já divergiram entre os ramos.

## Verificações numéricas e limites

- A suíte completa terminou com **438 testes aprovados e 3 ignorados**; não houve falhas. A suíte dirigida aos controles de evaporação, captura e análise teve 32/32 aprovações.
- O fechamento Euleriano máximo entre mudança de vorticidade e soma dos rotacionais dos incrementos foi de `5,95×10⁻¹⁸`, `6,75×10⁻¹⁸` e `6,15×10⁻¹⁸ s⁻²` em WEAK, CONTROL e STRONG.
- RK4 com passo máximo de 2 s versus 1 s diferiu por medianas de 0,0019–0,0022 m e máximos abaixo de 0,031 m. Usar saídas de 60 s em vez de 30 s mudou posições finitas por medianas de 6,32–6,39 m e máximos abaixo de 19,7 m.
- Foi imposta a sequência exata de `dt` do CONTROL para evitar que a lógica temporal virasse um segundo tratamento. Em WEAK, o `dt` imposto excedeu seu limite adaptativo próprio em no máximo 0,0924%; em STRONG, em no máximo 0,2304%. Nenhum passo excedeu seu limite próprio por 1%. Essa escolha preserva sincronização, mas é uma pequena assimetria numérica documentada.
- Existe uma única realização determinística por fator. O contraste controla exatamente as condições iniciais e o código, mas não estima dispersão estatística de um conjunto de perturbações.
- O primeiro evento de intensificação é censurado porque o reinício começa durante a intensificação. A perda final da conexão e o máximo completo de área fria são censurados pelo término em 3300 s.
- A malha horizontal de 600 m não resolve um tornado ou um funil estreito. O primeiro nível é um volume resolvido centrado acima do solo, não a superfície física. Todas as conclusões sobre rotação são da escala resolvida.
- A sequência compacta preserva o estado e os incrementos nos primeiros 2 km mais halo. A profundidade que atinge 2 km é apenas um limite inferior dentro dessa janela.
- Este teste não isola o cold pool como variável independente. Ele isola a evaporação da chuva; todas as respostas posteriores são mediação e feedback do modelo.

## Figuras e dados reproduzíveis

As 16 categorias obrigatórias e diagnósticos adicionais estão em `outputs/causal_evap_experiment/figures`:

1. `01_evaporacao_integrada.png`;
2. `02_area_cold_pool.png`;
3. `03_theta_v_minimo.png`;
4. `04_B_minimo.png`;
5. `05_zeta_baixa.png`;
6. `06_Tz.png`;
7. `07_stretching.png`;
8. `08_angulo.png`;
9. `09_orientacao_omega_h.png`;
10. `10_orientacao_grad_w.png`;
11. `11_componentes_grad_w.png`;
12. `12_mapas_horizontais_eventos.png`;
13. `13_cortes_verticais_evento_inversao.png`;
14. `14_trajetorias_lagrangianas.png`;
15. `15_comparacao_alinhada_por_eventos.png`;
16. `16_cadeia_causal.png`.

Os arquivos `17_profundidade_gradiente_posicao_cold_pool.png`, `18_efeitos_microfisicos_secundarios.png`, `19_velocidade_borda_equivalente.png` e `20_resposta_estado_strong_menos_weak.png` completam a auditoria. Séries numéricas estão em `final_metrics.csv`, `state_contrasts.csv`, `ordering_by_frame.csv`, `secondary_microphysics.csv` e nos CSVs de superfície/processo de cada ramo. Cada diretório de caso contém ainda `vortex_track.csv`, `parcels.csv`, `lagrangian_budget.csv`, `parcel_geometry.csv`, `orientation_budgets.csv` e metadados próprios.

## COMPROVADO PELO EXPERIMENTO

- O único tratamento deliberado foi `rain_evaporation_factor = 0,95/1,00/1,05`; CONTROL reproduziu o histórico bit a bit.
- O fator alterou de forma ordenada a evaporação efetivamente aplicada, o resfriamento latente, `qv`, `qr`, `θv`, `B`, a área e a intensidade média do cold pool.
- A inversão de `Tz` continuou sendo geometricamente dominada pela rotação de `∇h w`, e não pela rotação de `ωh`, em 27/27 parcelas de cada caso.
- O pico baixo, o pico material e o bracket de inversão caíram nas mesmas saídas nos três casos.
- A resposta de `ζ` após a inversão não seguiu H1/H2: STRONG manteve `ζ` ligeiramente maior e WEAK atingiu antes o critério de perda persistente.

## FORTEMENTE SUPORTADO

- A intervenção de ±5% é suficientemente pequena para preservar o regime geral e suficientemente grande para medir a primeira metade da cadeia causal.
- A evaporação da chuva causa uma resposta mensurável do cold pool nesta fase da simulação.
- O fator modula levemente `|∇h w|`, a rotação de `∇h w` e o mínimo de `Tz`, mas essa modulação não controla o instante resolvido da inversão.
- A hipótese de que maior evaporação, por meio de um cold pool mais intenso, antecipa a perda de `ζ` não é sustentada por este piloto.

## SUGESTIVO

- O cold pool mais intenso pode reforçar simultaneamente outras partes da circulação de baixo nível, pois STRONG apresentou `ζ` um pouco maior apesar de `Tz` mínimo um pouco mais negativo. Este teste não separa qual caminho produz essa compensação.
- A diferença de 0,76° na rotação de `∇h w` e de 1,28% no mínimo de `Tz` sugere uma sensibilidade geométrica real, mas pequena diante da evolução comum dominante.
- Os efeitos secundários microfísicos permanecem modestos e majoritariamente ordenados, compatíveis com feedback posterior à intervenção.

## NÃO DETERMINADO

- Um eventual deslocamento do horário da inversão menor que 30 s.
- O efeito de aplicar o fator desde o início da tempestade.
- O início original da intensificação, a perda final da conexão baixa e o máximo completo do cold pool.
- A profundidade do cold pool além de 2 km nas colunas que atingem o teto diagnóstico.
- O mecanismo que faz a pequena resposta de `ζ` ter sinal oposto ao esperado; não é atribuído aqui a LES, arrasto, baroclinicidade ou difusão numérica.
- Robustez estatística entre membros perturbados ou outras realizações.
- Formação de tornado ou de funil visível, que não é resolvida pela malha de 600 m e não foi usada como critério.

O experimento termina neste diagnóstico. Nenhum parâmetro adicional foi ajustado e nenhuma mudança física é proposta aqui.
