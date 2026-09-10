# Continuidade da análise de tornadogênese

## Reparação validada e proveniência em execução — 9 de setembro, atualização posterior

A simulação de 300 m terminou em 3300,023494218 s, passo 7397. O primeiro
pós-processamento falhou porque `DiagnosticCapture.save()` somava um ao
contador também quando chamado por `close()`, após o incremento do driver.
Somente o atributo `snapshots/00017.attrs.step` estava incorreto (7398).

Correção: `save(..., completed_step=None)` usa por padrão o contador atual;
o callback durante o passo passa explicitamente `sim.step+1`. `close()`
preserva os incrementos e escreve o contador já concluído. Teste explícito
CPU/GPU verifica os snapshots de passos 0, 3 e 4, sendo o último criado
apenas por `close()`. Resultado: 34 testes de captura/proveniência passaram.

`scripts/repair_capture_final_step.py` validou os 18 snapshots, formas,
valores finitos, índices 6114–7396, continuidade da agenda e soma dos dt:
1283 passos, 509,4439430736403 s, de 2790,579551144353 a 3300,023494218 s.
Reparou o atributo final para 7397 e registrou `metadata_corrections` no
HDF5. Hashes de todos os datasets antes/depois são idênticos; tempos,
agenda nativa e campos físicos não foram alterados.

SHA-256 do arquivo antes:
`40680e4995b55ba76d86c5fb819f5345768ea6007a4335eb357389dba9c4e736`.
Depois:
`2ade1b0911b453bba1fa808dc3871a8dc4e4f8bad612effa35b13c3368ec8fff`.
Registro completo: `outputs/resolution_300m_20260909/final_step_repair.json`.

Proveniência iniciada com os argumentos requeridos (snapshot 0 → 17, GPU,
`outputs/resolution_300m_provenance_20260909`). Não houve nova simulação
desde t=0. O runner agora registra a grade real e aceita diretório existente
somente quando vazio, para recuperar a tentativa anterior sem sobrescrever.
Limites originais de fechamento são verificados em cada estágio; neutralidade
bit a bit é verificada em cada passo. Aos 50/1283 passos: erro relativo de
ω 1,379e-15; κω=1,03327; κT=1,08378. **Gate final ainda pendente.**

`scripts/finish_repaired_resolution.py` está aguardando a proveniência;
verifica fechamento, finitude, neutralidade e restos advectivos, gera os
diagnósticos de 300 m com bins radiais comuns de 600 m e os artefatos de
comparação. Estado em `repaired_followup_status.json` no diretório da
sequência. Condicionamento e interpretação científica exigem revisão final.

Referência de 600 m reprocessada com as mesmas definições ampliadas em
`outputs/resolution_600m_comparison_audit_20260909`. Foram acrescentadas
integrais assinada/absoluta e máximo em cilindros de raio 4,2 km entre
0–2 km, circulação absoluta e larguras em metros. Comparação preparada em
`scripts/compare_resolution_audits.py`, com interpolação somente de métricas
escalares nos tempos comuns, sem extrapolação. Orçamentos legados em
máscaras conectadas móveis não equivalem a orçamento em volume fixo.

Auditoria geométrica de 300 m pronta em `outputs/domain_extent_300m_audit`:
contato de condensado com laterais e região superior também presente.
Permanece a classificação de sensibilidade conjunta à resolução e aos
efeitos dependentes do dt. Não foi demonstrada irrelevância do amortecimento.

**Resposta científica 300/600 m ainda não determinada:** concluir a execução,
revisar todos os gates e κ(t), conferir volumes/orçamentos e classificar
evidências antes de afirmar redução do déficit de concentração.

## Auditoria posterior do domínio (9 de setembro)

Ver `docs/DOMAIN_EXTENT_AUDIT.md`: nos 18 estados da janela madura de
600 m, condensado alcança células junto às laterais periódicas em todos
os estados e o último nível vertical em 16/18 para >0,01 g/kg. Cold pool
fica a pelo menos 11,1 km das bordas. Há suporte para investigar extensão
e topo, mas não prova causal de limitação da tornadogênese. O amortecimento
superior é aplicado por chamada, sem dt, e tem perfil invertido em relação
à rampa usual crescente até o teto; auditar antes de interpretar diferenças
entre resoluções como exclusivamente espaciais. Nenhuma física alterada.

## Execução refinada iniciada em 9 de setembro de 2026

O caso único de 300 m está em `outputs/resolution_300m_20260909`.
Inicialização analítica desde t=0, 240×240×48, mesma física e domínio.
Captura começa no primeiro passo nativo em ou após 2790,253490277 s;
fim em 3300,023494218 s. A diferença efetiva de tempo inicial deve ser
registrada na comparação. Não houve interpolação do estado de 600 m.
O spin-up não salva campos completos para reservar espaço em disco.

Comando: `python scripts/run_diagnostic_sequence.py --out outputs/resolution_300m_20260909 --nx 240 --nz 48 --duration 3300.023494218 --capture-start 2790.253490277 --interval 30 --device gpu`.

`scripts/finish_resolution_300m.py` foi iniciado como processo auxiliar:
aguarda `metadata.json` completo, executa proveniência v4 com controle
lado a lado e gera os diagnósticos somente após aprovação do gate.
Estado e erros ficam em `followup_status.json` e `followup.stderr.log`
dentro do diretório da execução. Saídas previstas:
`outputs/resolution_300m_provenance_20260909` e
`outputs/resolution_300m_audit_20260909`.

Validação antes do lançamento: 34 testes de captura/proveniência passaram.
O runner de proveniência agora lê nx/nz da sequência, mantendo a validação
exata da grade e do estado base no reinício.

**Pendente:** concluir a execução, verificar gates e comparar os produtos
de 300 e 600 m em volumes físicos comuns. Não há resultado científico
novo nem conclusão de convergência de malha nesta etapa. O auxiliar gera
os diagnósticos refinados; a comparação científica permanece pendente.

Atualizado em 8 de setembro de 2026. Este arquivo é autocontido para continuação em outra sessão.

## Estado científico atual

Classificação vigente: **MIXED — deficiência de concentração dominante, perda de inventário de circulação secundária**.

A tempestade mantém rotação forte aloft, produz vorticidade horizontal e executa tilting/stretching, mas não concentra uma coluna compacta e alinhada perto da superfície. Em 121,7 m, `ζmax` cai 29,7% do pico ao fim da janela e `Γ(4,2 km)` cai 56,0%. O máximo de convergência ocorre cerca de 240 s depois de `ζmax`; o eixo desloca-se 1,3–2,2 km entre 122 m e 2 km e até 5,16 km em níveis intermediários. O núcleo à meia amplitude tem 3,4–4,1 células na malha de 600 m. Isso não é um tornado resolvido.

## Métodos validados

1. `src/storm_dynamics/vorticity_provenance.py`, schema `storm-vorticity-provenance-v4`.
2. Parcelas de velocidade escalonada satisfazem `U=ΣUj`.
3. O operador de transporte usa a geometria C-grid MUSCL nativa com advector, sinal upwind e ramo minmod congelados a partir do total.
4. `ωj=curl(Uj)` é calculado depois de cada estágio. O `advection_remainder` guarda apenas arredondamento medido.
5. A opção `write_full_vorticity=True` salva `ξj`, `ηj` e `ζj`; é somente I/O e foi coberta pelo teste de fechamento.
6. O orçamento Euleriano primário usa o rotacional dos incrementos nativos capturados em todos os passos. A decomposição cinemática de transporte, stretching, tilting e dilatação é mantida separada.
7. O vórtice é rastreado em 3D a partir de 121,7 m por associação espacial contínua; máscaras conectadas usam `ζ≥0,002`, com robustez em `0,003` e `0,005 s⁻¹`.

Testes relevantes:

```powershell
python -m pytest tests/test_vorticity_provenance.py tests/test_vorticity_provenance_propagator.py -q
```

Resultado atual: `28 passed`.

## Métodos inválidos ou retraídos

- O propagador antigo de ω com Euler explícito e diferenças centradas é inválido. Em cerca de 2981 s atingiu `κω≈3,98×10⁴`, `κT≈2,29×10⁵`.
- Todos os produtos físicos anteriores marcados por `outputs/vorticity_provenance_analysis/INVALID_FOR_CURRENT_PROTOCOL.md` permanecem retraídos e não foram usados.
- Não existe trajetória Lagrangiana validada contra o transporte MUSCL. Não reutilizar o integrador Euler centrado antigo e não alegar proveniência parcelar.
- Não reinterpretar v4 como Jacobiano completo. É uma convenção aditiva frozen-advector/frozen-limiter.
- Não reutilizar a interpretação Nyquist retraída de 26 m/s.

## Execução decisiva

- sequência total: `outputs/diagnostic_sequence_20260905/sequence.h5`;
- reinício: snapshot `00093`, `t=2790,253490277 s`, passo `4406`;
- alvo: snapshot `00110`, `t=3300,023494218 s`, passo `5421`;
- intervalo: `509,770004 s`, 1.015 passos nativos;
- grade: `120×120×48`, `Δx=Δy=600 m`, vertical esticada;
- dispositivo: NVIDIA GeForce RTX 4050 Laptop GPU;
- fator de evaporação: `1,0`;
- commit de base: `097d516409285efad047012f80f8b2b8c615f71b`;
- working tree: modificada; hashes exatos das fontes estão no `summary.json` da execução.

Comando:

```powershell
python scripts/run_vorticity_provenance_long.py --snapshot 93 --target-snapshot 110 --device gpu --interval 30 --out outputs/vorticity_provenance_long_v4_2790_3300
```

Resultado do gate:

- todos os 12 prognósticos bit a bit idênticos;
- `max RMS(Rω)/RMS(ω)=5,78×10⁻¹⁵`;
- `max|Rω|=7,53×10⁻¹⁶ s⁻¹`;
- `max RMS(RT)/RMS(T)=6,89×10⁻¹⁵`;
- `κω(final)=1,851`, `κT(final)=2,691`;
- restos advectivos relativos finais: `6,55×10⁻¹⁵` para ω e `7,01×10⁻¹⁵` para tilting;
- HDF5: `1.592.700.738 bytes`, SHA-256 `a115ed2b1dae3497fa7757ce89e1744db83f5ad0139e8c10989a9b254e5a3ccf`;
- tempo de parede: `1.981,4 s`.

O atributo `step` do quadro final foi corrigido de 5422 para 5421 depois da execução. O próprio HDF5 registra essa correção em `metadata_corrections`; nenhum array ou tempo físico mudou. A lógica de `close()` foi corrigida e ganhou teste de regressão.

## Análise reproduzível

Comando:

```powershell
python scripts/analyze_tornadogenesis_mechanism.py
```

Saídas:

- `outputs/tornadogenesis_mechanism_audit/summary.json`;
- `manifest.json`, com metadados e SHA-256 dos produtos;
- `vortex_timeseries.csv`;
- `threshold_sensitivity.csv`;
- `vertical_profiles.csv`;
- `radial_profiles.csv`;
- `provenance_integrals.csv`;
- `vorticity_budget.csv`;
- `principal_events.csv`;
- `event_states.h5`, com os campos compactos 0–2 km de cada tempo de evento;
- `lag_correlations.csv`;
- `figures/`, com 32 PNGs.

A correspondência entre quadros totais e de proveniência difere no máximo `0,576 s`. O fechamento de proveniência é calculado nos tempos exatos; essa correspondência só afeta sobreposições com termodinâmica e pressão.

## Resultados numéricos principais

### Vórtice de baixo nível

- `ζmax(121,7 m)`: início `0,00533`, pico `0,00554` em 2880 s, final `0,00390 s⁻¹`.
- `Γ(4,2 km)`: início/pico `34.437`, final `15.160 m² s⁻¹`.
- `Vθ,max`: máximo `2,28 m s⁻¹`.
- queda de `p_dyn`: mínimo `−45,85 Pa` em 3151 s.
- convergência máxima: `0,01040 s⁻¹` em 3120 s.
- `wmax` em 121,7 m: máximo `1,315 m s⁻¹` em 3060 s.
- RMW: 0,5–4,5 células; abaixo de duas células em 10/18 quadros.
- largura de ζ à meia amplitude: 3,4–4,1 células.
- largura de convergência à meia amplitude: 8,1–9,8 células.
- largura da pressão à meia queda: 7,5–7,9 células.

### Proveniência em 3300 s

- ζ assinado: `initial 63,42%`, `projection 36,92%`, `coriolis +2,39%`, `les −2,70%`, drag `−0,032%`.
- norma horizontal: `initial 42,72%`, `projection 28,70%`, `buoyancy 20,27%`, `les 4,56%`.
- tilting assinado: `projection 78,81%`, `initial 14,62%`, `buoyancy 13,01%`, `coriolis −5,90%`, `les −2,13%`.
- `advection_remainder` é numericamente desprezível.

### Orçamento Euleriano

Médias temporais na máscara móvel: mudança total `−4,77×10⁻⁷ s⁻²`, advecção MUSCL `−2,04×10⁻⁷`, LES `−5,36×10⁻⁷`, arrasto `−2,54×10⁻⁸`, Coriolis `+2,89×10⁻⁷`, rotacional vertical direto da projeção `≈0`. Cinemática independente: stretching `+2,10×10⁻⁵`, tilting `+1,79×10⁻⁵`, transporte `−3,13×10⁻⁵`, dilatação `−4,60×10⁻⁶ s⁻²`.

### Cold pool

- área com `θv'<−1 K`: 169,9–290,5 km²;
- 37,4–44,2% da vizinhança de 4,2 km está fria;
- centro fica uma célula dentro ou fora da interface (`±600 m`);
- geração baroclínica horizontal média local: `1,54–1,76×10⁻⁵ s⁻²`.

O contato e a geração estão demonstrados. O sinal causal do cold pool permanece indeterminado.

## Interpretação e hipóteses remanescentes

1. **Concentração/alinhamento:** hipótese dominante, fortemente suportada por dessincronização, baixa correlação espacial, eixo inclinado e núcleo marginalmente resolvido.
2. **Perda de fonte/inventário:** componente secundário, comprovado pela queda de circulação e pelas integrais em todos os limiares.
3. **Conversão insuficiente:** não suportada como dominante; horizontal ω, tilting e stretching são abundantes.
4. **LES/difusão:** sumidouro associado quantitativamente, mas causalidade não determinada.
5. **Cold pool:** fonte e interface presentes; intensidade ótima ou excessiva não determinada.
6. **História pré-2790:** não separada, pois está em `initial`.

## Próximo passo exato

Executar **um único caso de resolução `Δx=300 m`**, com estados iniciais/reinício dinamicamente compatíveis e toda a física idêntica. Usar v4 na mesma janela e comparar `Γ(r)`, `Vθ(r)`, RMW em células, larguras de ζ/convergência/pressão, alinhamento vertical e orçamento LES/MUSCL em volumes físicos comuns. Não iniciar varredura de LES, drag, cold pool ou evaporação antes desse teste.

Relatório completo: `docs/TORNADOGENESIS_MECHANISM_AUDIT.md`.

## Atualização 2026-09-09 — caso de 300 m concluído

A proveniência v4 percorreu os 1.283 passos nativos entre 2790,579551 e
3300,023494 s. O gate numérico passou, os 12 campos prognósticos permaneceram
bit a bit idênticos e o máximo fechamento relativo foi `8,7053e-15`.

A comparação final usa 18 tempos comuns, nível de 121,659 m, bins radiais de
600 m, cilindro de raio 4,2 km e volume entre 0 e 2 km. As quantidades em células
foram separadas das mudanças físicas após revisão independente.

Em unidades físicas, 300 m elevou a mediana de zeta máxima de baixo nível em
3,10 vezes, a circulação de 4,2 km em 32% e a velocidade tangencial em 31%; a
largura de zeta à meia altura caiu 46%. O déficit de concentração diminuiu, mas
o RMW ainda possui cerca de três células e o vórtice continua enfraquecendo. A
classificação permanece `MIXED: CONCENTRATION DEFICIT DOMINANT,
SOURCE-INVENTORY LOSS SECONDARY`.

O domínio de 72 x 72 x 15 km mantém o cold pool a mais de 11 km das bordas, mas
o condensado alcança a costura periódica e o topo. O damping superior é aplicado
quatro vezes por passo, sem escala por dt, então a experiência é uma sensibilidade
conjunta à resolução e a efeitos numéricos dependentes do passo de tempo.

Relatório completo: `docs/RESOLUTION_300M_COMPARISON.md`.

## Atualização 2026-09-10 — extensão lateral para 120 km concluída

O experimento isolado ampliou apenas `Lx`, `Ly`, `nx` e `ny`, de
72×72×15 km (`120×120×48`) para 120×120×15 km (`200×200×48`), mantendo
`dx=dy=600 m`, a malha vertical, a física e o movimento da tempestade. O
preflight confirmou identidade bit a bit dos 12 campos prognósticos na região
central comum e um passo real finito na GPU.

A simulação terminou em `3300,023494218 s`, com 5.578 passos, 18 quadros e
1.103 passos nativos na janela capturada. O arquivo `sequence.h5` tem SHA-256
`21169b14beb545d4e21698fe73b8d195e65b05272e71d4dfc4fa039584ca4513`.

O gate lateral passou. A distância mínima à borda foi 6,3 km para condensado
leve (`qcond>1e-5`) e 26,7 km para corrente ascendente acima de 5 m/s. Portanto,
120 km são suficientes para conter horizontalmente a supercélula durante a
janela analisada, embora o limiar de condensado leve passe com margem de apenas
0,3 km.

A ampliação lateral não melhorou materialmente o vórtice: 0/5 métricas
primárias cumpriram simultaneamente coerência em pelo menos 12/18 pares e efeito
mediano de pelo menos 20%. Em relação a 72 km, as medianas mudaram em `-3,89%`
para zeta máxima de baixo nível, `-5,89%` para circulação a 4,2 km, `-3,46%`
para velocidade tangencial máxima, `+4,53%` para largura de zeta e `-2,56%`
para inclinação do eixo entre a superfície e 2 km. A classificação
pré-registrada é `no dominant favorable domain effect`.

O gate vertical permanece `UNRESOLVED`: condensado leve alcança o centro da
última célula (`14,605 km`), com até `41,04 km²` no último nível, e há
movimento vertical de `1,05–1,34 m/s` na região amortecida. Além disso, o caso
de 120 km recebeu 8,65% mais exposição nominal ao damping superior por segundo
na janela, devido ao passo adaptativo; esse é um confundidor de realimentação.

Conclusão atual: o domínio original de 72 km era pequeno para a nuvem completa,
mas isso não explica o déficit de intensidade do vórtice observado. O próximo
teste geométrico adequado é estender o topo para 20 km preservando exatamente
as faces da malha abaixo de 15 km; mudar `Lz` mantendo `nz=48` confundiria altura
e resolução vertical. Uma afirmação causal forte ainda requer ensemble pareado.

Relatório completo: `docs/LATERAL_DOMAIN_RESULTS.md`.
