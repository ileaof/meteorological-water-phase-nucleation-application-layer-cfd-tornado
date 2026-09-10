# Auditoria do mecanismo de tornadogênese

Data: 8 de setembro de 2026.

## 1. Conclusão executiva

**Classificação: MIXED, com CONCENTRATION DEFICIT dominante e perda de inventário associada a SOURCE DEFICIT secundário.**

O vórtice de baixo nível já existe no início da janela, mas não contrai nem mantém sua circulação. Em 121,7 m, `ζmax` atinge `0,00554 s⁻¹` em 2880 s e cai 29,7% até 3300 s; a circulação em 4,2 km cai 56,0%, de `34.437` para `15.160 m² s⁻¹`. Há vorticidade horizontal e tilting abundantes, portanto os dados não sustentam uma deficiência primária de conversão. O stretching é positivo em 82,5% das células de vorticidade positiva, mas a correlação mediana entre ζ e convergência é `−0,045`, o máximo de convergência ocorre cerca de 240 s depois de `ζmax`, e o eixo sofre deslocamentos de 1,3–2,2 km entre 122 m e 2 km. O núcleo de ζ tem somente 3,4–4,1 células à meia amplitude, enquanto convergência e pressão ocupam 8–10 células. O resultado é um mesociclone de baixo nível amplo, inclinado e intermitentemente conectado, sem concentração resolvida, vento tangencial ou queda dinâmica de pressão compatíveis com um vórtice tornádico inequívoco.

## 2. Validade numérica

A execução autorizada usou o propagador `storm-vorticity-provenance-v4`, de 2790,253490 a 3300,023494 s, em 1.015 passos nativos na GPU, com uma solução sem instrumento avançada lado a lado.

### Neutralidade

`u`, `v`, `w`, `p`, `theta`, `qv`, `ql`, `qi`, `qr`, `qs`, `qg` e `qh` permaneceram **idênticos bit a bit após todos os passos**. Para cada variável:

- diferença absoluta máxima: `0`;
- diferença RMS máxima: `0`;
- diferença RMS relativa máxima: `0`.

### Fechamento e condicionamento

| diagnóstico, 0–2 km | máximo na janela | final |
|---|---:|---:|
| `RMS(Rω)/RMS(ω)` | `5,78×10⁻¹⁵` | `5,77×10⁻¹⁵` |
| `max|Rω|` | `7,53×10⁻¹⁶ s⁻¹` | `6,71×10⁻¹⁶ s⁻¹` |
| `RMS(RT)/RMS(T)` | `6,89×10⁻¹⁵` | `6,86×10⁻¹⁵` |
| `max|RT|` | `2,53×10⁻¹⁸ s⁻²` | `2,18×10⁻¹⁸ s⁻²` |
| `||Radv,ω||/||ω||` | `6,55×10⁻¹⁵` | `6,55×10⁻¹⁵` |
| `||Radv,ω||/Σ||ωj||` | `3,54×10⁻¹⁵` | `3,54×10⁻¹⁵` |
| `||Radv,T||/||T||` | `7,01×10⁻¹⁵` | `7,01×10⁻¹⁵` |
| `||Radv,T||/Σ||Tj||` | `2,60×10⁻¹⁵` | `2,60×10⁻¹⁵` |
| `κω` | `1,851` | `1,851` |
| `κT` | `2,696` | `2,691` |

Os índices crescem suavemente a partir de 1, sem crescimento secular ou exponencial. Permanecem muito próximos da validação madura curta (`κω=1,267`, `κT=1,655` em 150 s) e várias ordens de grandeza abaixo do método rejeitado (`3,98×10⁴` e `2,29×10⁵`). O orçamento formado pelo rotacional de todos os incrementos nativos fecha a mudança de ζ com erro RMS relativo máximo de `3,29×10⁻¹⁵`.

O gate numérico foi **APROVADO**. A interpretação abaixo não usa os produtos físicos v1/v2/v3 retraídos.

## 3. Evolução do vórtice de baixo nível

O objeto foi rastreado em 3D a partir do máximo ciclônico em 121,7 m. Em cada nível, o centro é a média ponderada por ζ positiva ao redor do máximo associado ao centro anterior; a associação vertical admite deslocamento máximo de 4,2 km por nível. O volume principal é o componente conectado a esse eixo com `ζ≥0,002 s⁻¹`. A robustez foi verificada com `0,003` e `0,005 s⁻¹`, além das métricas contínuas de circulação, pressão, vento tangencial, convergência e `w`.

- **2790 s:** o vórtice já está presente; a janela não contém o início de sua formação. A circulação em 4,2 km já está em seu maior valor observado, `34.437 m² s⁻¹`.
- **2880 s:** `ζmax` em 121,7 m chega a `0,00554 s⁻¹`, apenas 4,0% acima do início.
- **2941 s:** `Vθ,max` chega a `2,28 m s⁻¹`, mas o RMW estimado é 300 m, meia célula, portanto esse raio não é resolvido. A circulação já caiu para `25.949 m² s⁻¹`.
- **3060 s:** ocorre o enfraquecimento mais rápido de `ζmax`, apesar de `wmax=1,31 m s⁻¹` e convergência local forte.
- **3120 s:** a convergência máxima, `0,01040 s⁻¹`, ocorre depois da intensificação de ζ.
- **3151 s:** a maior queda de pressão dinâmica local, `−45,85 Pa`, ocorre quando a circulação caiu para `15.900 m² s⁻¹`.
- **3300 s:** `ζmax=0,00390 s⁻¹`, circulação `15.160 m² s⁻¹`, `Vθ,max=1,97 m s⁻¹` e queda dinâmica `−25,94 Pa`.

O evento “início do crescimento de ζ” é censurado pelo limite esquerdo: o vórtice já existe em 2790 s. As correlações defasadas usam somente 18 quadros e são apresentadas como exploratórias, não como causalidade.

## 4. Proveniência da vorticidade

A decomposição é a convenção aditiva validada, com advector, sinais upwind e ramos do limiter congelados no estado total. `initial` significa toda a história anterior ao reinício; não é um mecanismo físico único.

No volume rastreado em 3300 s, as frações do **inventário assinado de ζ** são:

| rótulo | fração assinada |
|---|---:|
| `initial` | `63,42%` |
| `projection` | `36,92%` |
| `coriolis` | `+2,39%` |
| `les` | `−2,70%` |
| `surface_drag` | `−0,032%` |
| `buoyancy`, `boundary`, `other` | `0%` direto |
| `advection_remainder` | `<4×10⁻¹⁵%` |

A contribuição direta de `buoyancy` para ζ é zero porque a força é aplicada a `w`; seu rotacional gera inicialmente vorticidade horizontal. Isso não significa ausência de participação na conversão. No inventário de `|ωh,j|` em 3300 s, `initial` representa 42,72%, `projection` 28,70%, `buoyancy` 20,27%, `les` 4,56%, `surface_drag` 1,57%, `coriolis` 1,42% e `boundary` 0,75%.

No tilting condicionado e integrado, as frações assinadas finais são `projection=78,81%`, `initial=14,62%`, `buoyancy=13,01%`, `coriolis=−5,90%`, `les=−2,13%`, `surface_drag=+1,23%` e `boundary=+0,36%`. `projection` tem rotacional vertical praticamente nulo no instante em que é aplicado, mas as velocidades horizontais que ele cria são depois transportadas e deformadas. Esse é um exemplo concreto da diferença entre proveniência acumulada e termo Euleriano instantâneo.

O fator `initial` ainda responde por 59,78% do inventário positivo normalizado em 3300 s. Assim, esta janela prova a evolução pós-2790, mas não atribui retroativamente a origem física de toda a vorticidade antecedente.

## 5. Orçamento de vorticidade

O orçamento primário usa o rotacional das mudanças de velocidade realmente aplicadas em cada estágio do solver, acumuladas em cada intervalo de aproximadamente 30 s. As médias temporais no volume rastreado são:

| termo nativo | média (`s⁻²`) |
|---|---:|
| mudança total de ζ | `−4,77×10⁻⁷` |
| advecção MUSCL líquida | `−2,04×10⁻⁷` |
| LES | `−5,36×10⁻⁷` |
| arrasto superficial | `−2,54×10⁻⁸` |
| Coriolis | `+2,89×10⁻⁷` |
| projeção, rotacional vertical direto | `−3,1×10⁻²¹` |
| demais estágios | zero na precisão armazenada |

O sinal líquido negativo é reproduzido pelo orçamento discreto. LES é o maior termo negativo médio dentro da máscara móvel, mas este diagnóstico não demonstra que LES cause a falha: mudança de máscara, resolução e transporte estão acoplados, e não foi feito um experimento de sensibilidade.

A decomposição cinemática independente fornece stretching médio `+2,10×10⁻⁵ s⁻²`, tilting `+1,79×10⁻⁵ s⁻²`, transporte de vorticidade `−3,13×10⁻⁵ s⁻²` e dilatação `−4,60×10⁻⁶ s⁻²`. Esses termos contínuos diagnosticados não devem ser somados aos termos nativos como fontes adicionais; eles decompõem geometricamente a advecção/deformação cuja atualização efetiva é o termo MUSCL. Mostram que produção local positiva é largamente compensada por remoção e redistribuição.

Não existe no solver um termo vertical baroclínico aplicado separadamente. A flutuabilidade vertical gera vorticidade horizontal por gradientes horizontais de densidade, posteriormente disponível para tilting.

## 6. Stretching e convergência

Não falta stretching em magnitude:

- média condicional em `ζ≥0,005 s⁻¹`: `2,60–3,03×10⁻⁵ s⁻²`;
- mediana da fração de células com ζ positiva e stretching positivo: `82,51%`;
- mediana da fração com ζ positiva e convergência positiva: `81,58%`.

O problema é a organização desse forcing. `corr(ζ,C)` tem mediana `−0,045` e varia somente de `−0,115` a `+0,072`; `corr(ζ,∂w/∂z)` tem mediana `+0,059`. A convergência máxima permanece perto de `0,010 s⁻¹`, mas o seu núcleo à meia amplitude tem 8,1–9,8 células, enquanto o núcleo de ζ tem 3,4–4,1. A queda de pressão é igualmente ampla, 7,5–7,9 células. Convergência e stretching atuam em grande parte do objeto, porém não se concentram e sincronizam no núcleo ciclônico de baixo nível.

## 7. Circulação versus `ζmax`

O resultado rejeita a interpretação de “contração com circulação conservada”:

- `ζmax`: `0,00533 → 0,00554 → 0,00390 s⁻¹`, final 70,3% do pico;
- `Γ(4,2 km)`: `34.437 → 15.160 m² s⁻¹`, final 44,0% do pico;
- raio equivalente do componente `ζ≥0,002 s⁻¹`: aproximadamente `1,35 → 1,02 km`.

O objeto encolhe, mas ζ não aumenta; ao contrário, o inventário de alta ζ diminui. A integral condicionada cai 8,35% para o limiar `0,002`, 11,95% para `0,003` e 27,82% para `0,005 s⁻¹`. Portanto há perda/redistribuição de circulação e falha de concentração simultâneas.

## 8. Coerência vertical

O ramo associado ao vórtice de baixo nível é contínuo pelo critério local, mas não é aproximadamente vertical:

- deslocamento entre 121,7 m e 1.974 m: `1,31–2,20 km`;
- maior deslocamento em qualquer nível de 0–2 km: `2,52–5,16 km`;
- em 2790 s, a circulação passa de `34.437 m² s⁻¹` em 122 m para `151.538 m² s⁻¹` perto de 2 km;
- em 3300 s, passa de `15.160` para `158.742 m² s⁻¹`.

A rotação aloft permanece muito mais forte que a rotação de baixo nível. Há um mesociclone verticalmente conectado por um caminho inclinado, mas a coluna de forte ζ não desce como um tubo compacto e alinhado. Uma falha de alinhamento faz parte da deficiência de concentração.

## 9. Interação com o cold pool

Usando `θv'=θv−θv0` e a borda `θv'=−1 K` em 121,7 m:

- área fria: `169,9–290,5 km²`;
- fração fria na vizinhança de 4,2 km do vórtice: `37,4–44,2%`;
- distância assinada do centro à borda: `−600` ou `+600 m`, isto é, uma célula dentro ou fora;
- geração baroclínica horizontal média perto do vórtice: `1,54–1,76×10⁻⁵ s⁻²`;
- `buoyancy` responde por 20,27% da soma das normas de vorticidade horizontal e 13,01% do tilting assinado em 3300 s.

Está comprovado que o vórtice acompanha a interface do cold pool e que esta fornece vorticidade horizontal utilizável. Não está determinado se o cold pool é excessivamente forte, fraco ou mal posicionado. Temperatura, proximidade e correlação não bastam para atribuir causalidade.

## 10. Resolução e difusão numérica

A malha efetivamente usada tem `Δx=Δy=600 m`, e não 444 m. O diâmetro do núcleo de ζ à meia amplitude ocupa somente 3,4–4,1 células. O RMW diagnosticado varia de 0,5 a 4,5 células e fica abaixo de duas células em 10 dos 18 quadros. Esses valores abaixo de duas células indicam falha do estimador em resolver o raio, não contração física demonstrada.

A convergência e a pressão são resolvidas em estruturas mais largas, mas a escala que precisaria contrair para um vórtice tornádico está no limite da malha. O termo LES é um sumidouro persistente no orçamento local (`−5,36×10⁻⁷ s⁻²` em média). Isso torna resolução/difusão um mecanismo limitante plausível e quantitativamente motivado. Ainda não permite separar LES explícita de dissipação MUSCL nem declarar causalidade sem uma comparação refinada.

## 11. Diagnóstico da tornadogênese

**Por que a tempestade simulada ainda não produz um vórtice tornádico inequívoco?**

Porque a cadeia alcança geração de vorticidade horizontal, tilting e stretching, mas não organiza esses ingredientes em uma concentração compacta, alinhada e sustentada de circulação perto da superfície. O pico de convergência chega depois do pico de ζ, os campos de convergência/pressão são largos em relação ao núcleo, o eixo se desloca por vários quilômetros ao longo dos primeiros 2 km e a circulação de baixo nível já está diminuindo. A resposta dinâmica resolvida permanece fraca: `Vθ,max≤2,28 m s⁻¹` e `|Δpdyn|≤45,9 Pa`. Esses números descrevem o vórtice resolvido; não são limiares observacionais de tornado.

Avaliação das hipóteses:

- **H1 — SOURCE DEFICIT: suportado como componente secundário.** O inventário/circulação de baixo nível diminui e nenhuma fonte pós-reinício reverte essa perda, embora haja vorticidade horizontal abundante.
- **H2 — CONVERSION DEFICIT: não suportado como limitação dominante.** Tilting e stretching são positivos e crescem ou permanecem fortes enquanto ζ enfraquece.
- **H3 — CONCENTRATION DEFICIT: fortemente suportado e dominante.** Há dessincronização, fraca correlação espacial, eixo inclinado, núcleo marginalmente resolvido e ausência de crescimento de ζ durante a contração aparente.

## 12. Incertezas remanescentes

1. A origem anterior a 2790 s permanece agregada em `initial` (63,42% do ζ assinado final).
2. Não há método de trajetórias numericamente validado contra o transporte MUSCL; proveniência material/parcelar permanece não determinada e o integrador Euler centrado antigo não foi reutilizado.
3. A convenção frozen-advector é aditiva e fechada, mas não é a decomposição física única da não linearidade.
4. A interpolação temporal entre sequência total e quadros de proveniência é no máximo 0,576 s; os fechamentos internos usam tempos exatos, enquanto sobreposições com termodinâmica usam o quadro total mais próximo.
5. O sinal do efeito causal de cold pool, LES, arrasto e difusão numérica não foi determinado. O orçamento quantifica associações, não substitui perturbações controladas.
6. O RMW e qualquer estrutura menor que aproximadamente 2–4 células não são quantitativamente convergidos.
7. `p_dyn` é a pressão de projeção diagnóstica e não é acoplada à saturação; queda dinâmica e funil visível permanecem diagnósticos separados.

## 13. Próximo experimento

Executar **uma única repetição refinada para `Δx=300 m`**, mantendo física, termodinâmica, hodógrafa, arrasto, LES, fator de evaporação (`1,0`) e intervalo de análise. A comparação deve começar de estados dinamicamente compatíveis, evitar uma interpolação de reinício não validada e usar o mesmo propagador v4 e as mesmas máscaras físicas. Os critérios decisivos são:

1. convergência de `Γ(r)` e `Vθ(r)` em raios físicos comuns;
2. aumento de `RMW/Δx` e largura do núcleo acima da faixa marginal atual;
3. redução ou persistência da perda de circulação;
4. mudança ou persistência da separação entre ζ, convergência, stretching e pressão;
5. orçamento LES/MUSCL por volume físico comum.

Se o refinamento produzir contração e pressão mais fortes com circulação comparável, a limitação numérica de concentração será confirmada. Se o padrão persistir em células e metros, a deficiência dinâmica de organização ganhará suporte. Não foi iniciada uma varredura de parâmetros.

## Artefatos

- execução v4 e gate: `outputs/vorticity_provenance_long_v4_2790_3300/summary.json`;
- campos 3D de proveniência: `outputs/vorticity_provenance_long_v4_2790_3300/provenance_long_v4.h5`;
- síntese: `outputs/tornadogenesis_mechanism_audit/summary.json`;
- metadados e hashes dos produtos: `outputs/tornadogenesis_mechanism_audit/manifest.json`;
- tabela de eventos: `outputs/tornadogenesis_mechanism_audit/principal_events.csv`;
- estados compactos 0–2 km dos eventos: `outputs/tornadogenesis_mechanism_audit/event_states.h5`;
- séries e orçamentos: `outputs/tornadogenesis_mechanism_audit/*.csv`;
- 32 figuras: `outputs/tornadogenesis_mechanism_audit/figures`.
