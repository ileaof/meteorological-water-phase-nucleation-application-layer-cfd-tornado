# Controle experimental da evaporação da chuva

**Existe agora um parâmetro contínuo suficientemente isolado para realizar o teste causal do cold pool? SIM, COM LIMITAÇÕES.**

Foi introduzido `rain_evaporation_factor`, adimensional, padrão **1.0**, exclusivamente na taxa de evaporação da chuva antes dos limitadores existentes. O valor padrão reproduz a implementação anterior bit a bit nos testes CPU/GPU e na amostra da sequência salva. Não foi executado nenhum caso WEAK-CP, CONTROL ou STRONG-CP, nem selecionada a amplitude de uma futura intervenção.

O isolamento é **da taxa microfísica**, não do cold pool como variável independente. A mesma transferência física muda chuva, vapor, calor latente e carga de condensado. Esses efeitos são acoplados pela conservação; não é possível mudar somente a temperatura sem modificar a física. O fator permite um experimento causal sobre a evaporação de chuva, com análise posterior da mediação pelo cold pool.

## Rotina e cadeia de chamadas

O ponto de intervenção é [`rain_evaporation`](../src/precip_microphysics/processes.py), no módulo `precip_microphysics.processes`. No código de produção, a chamada ocorre indiretamente pelo último elemento de `PROCESS_ORDER`, iterado por `BulkMicrophysics.step` em [`scheme.py`](../src/precip_microphysics/scheme.py). O esquema é chamado por `MicrophysicsCoupler.apply` na circulação 3D e por `ColumnModel` no modelo de coluna. `StormSimulation._transport` chama o acoplador; `meteorological_flow.Simulation` usa o mesmo caminho microfísico. Os novos testes e o script local chamam a rotina diretamente.

`rain_evaporation` não é chamada por gelo, neve, graupel ou hail. O que é compartilhado é `_ventilated_capacitance`, usada pela evaporação de chuva, pela deposição/sublimação de neve e por `_melt`, que atende derretimento de neve/graupel/hail. `_diffusional_denominator` também é compartilhada por processos difusionais. **Nenhuma dessas funções compartilhadas foi modificada.**

A evaporação de água de nuvem pertence a `condensation_adjustment`, cujo interruptor também controla condensação. Ela não recebe o novo fator. Autoconversão, acreção, congelamento, agregação, riming, sedimentação e nucleação tampouco o leem.

## Equação efetivamente implementada antes da alteração

A rotina retorna uma transferência de massa `Transfer("qr", "qv", dq, "rain_evaporation")`, com `dq` em kg/kg. Para o estado de entrada e passo Δt:

\[
\lambda=\left(\frac{\pi\rho_wN_{0r}}{\rho q_r}\right)^{1/4},\qquad
S_w=\frac{p_v(q_v,P)}{e_{sw}(T)}.
\]

`lambda_slope` usa a distribuição exponencial de momento único, com `N0_r=8.0e6 m⁻⁴`, `rho_w=1000 kg/m³`. Os campos opcionais de número de gotas `Nr` não entram nessa taxa: o tamanho é diagnosticado de qr, ρ e N0. A pressão de vapor usa a convenção de fração mássica do repositório; e_sw é calculada por `thermo.psat_water`.

\[
A_K=\frac{L_v}{K_{THERM}T}\left(\frac{L_v}{R_vT}-1\right),\qquad
A_D=\frac{R_vT}{D_v\max(e_{sw},TINY)}.
\]

O argumento P de `_diffusional_denominator` não é usado dentro dessa função. P influencia a taxa por S_w e o limitador q_sat; T influencia S_w, e_sw e a resistência difusional. A densidade entra na distribuição e no prefator. Não foi adicionada uma dependência de pressão ao coeficiente de difusão.

A capacitância ventilada integrada é calculada exatamente como:

\[
V=N_{0r}\left[\frac{0.78}{\lambda^2}+
0.31\,Sc^{1/3}\sqrt{a/\nu_{air}}\,
\Gamma\!\left(\frac{b+5}{2}\right)\lambda^{-(b+5)/2}\right],
\]

com `a=841.99667`, `b=0.8` para chuva, `NU_AIR=1.5e-5 m²/s`, `DIFF_VAPOR=2.26e-5 m²/s`, `Sc=NU_AIR/DIFF_VAPOR` e `K_THERM=2.43e-2 W m⁻¹ K⁻¹`. As unidades de a seguem `Vt=a D^b`, portanto dependem de b.

\[
r=\frac{2\pi}{\max(\rho,TINY)}\frac{S_w-1}{A_K+A_D}V,
\qquad R_0=\begin{cases}-r,&r\text{ finito e }S_w<1\\0,&\text{caso contrário.}\end{cases}
\]

R₀ é a magnitude não negativa de evaporação, em kg kg⁻¹ s⁻¹. Antes disso, a rotina retorna sem transferência quando o interruptor está desligado ou não há nenhuma célula subsaturada com `qr > QSMALL`. `QSMALL=1e-12 kg/kg`; a distribuição também mascara categorias abaixo desse limiar.

Finalmente:

\[
\Delta q_r^{evap}=\operatorname{clip}
\left[\min(R_0\Delta t,\max(q_{sat,w}(T,P)-q_v,0)),\;0,\;\max(q_r,0)\right].
\]

O limitador de saturação usa T de entrada; a transferência aquece/resfria o estado posteriormente no esquema. Essa ordem e os limitadores foram preservados, sem corrigir ou reinterpretar sua física.

### Constantes declaradas versus expressão usada

`constants.VENT_A=0.78` e `VENT_B=0.31` existem e são usados por `thermo.ventilation_factor`. Entretanto, o caminho de taxas acima não chama essa função: `_ventilated_capacitance` contém os literais 0.78 e 0.31. Portanto, mudar somente VENT_A/VENT_B não controla esta evaporação. Eles não são universalmente “não utilizados”; estão desconectados especificamente deste caminho.

Os literais, as constantes físicas, a distribuição de gotas, o cálculo de saturação e todos os coeficientes foram mantidos. Não houve deslocamento de coeficientes para outro lugar nem refatoração numérica da ventilação.

## Alteração mínima e configuração

A única mudança da equação de intervenção é:

\[
R_{new}=fR_0,\qquad f=\texttt{rain_evaporation_factor}.
\]

O fator é aplicado a `evap`, imediatamente antes de `min(evap * dt, to_sat)`. Para f=1, o ramo de multiplicação é pulado e a sequência aritmética anterior permanece intacta. Para f=0, a rotina de chuva retorna sem transferência. O interruptor antigo continua funcionando; f=0 não desliga evaporação de nuvem nem sublimação de gelo.

`MicrophysicsConfig` valida um escalar finito e não negativo; valores negativos, NaN e infinito são rejeitados. O mesmo valor é exposto em `SimulationConfig.physics`, encaminhado pelos dois drivers ao acoplador e incluído nos metadados resumidos do driver de circulação. A serialização completa por `asdict`, usada pela instrumentação da tempestade, também o inclui.

Configuração microfísica independente:

```yaml
rain_evaporation_factor: 1.0
```

Configuração do escoamento:

```yaml
physics:
  rain_evaporation_factor: 1.0
```

No objeto de tempestade, o caminho é `cfg.sim.physics.rain_evaporation_factor`, definido **antes** de construir `StormSimulation`. Arquivos antigos sem a chave assumem 1.0 automaticamente. O fator é copiado ao construir o acoplador; não se deve alterar apenas o objeto externo depois da construção e esperar atualização automática.

Não foi criado um ramo diferente para CPU, CuPy ou projeção de baixa memória. O escalar é aplicado à matriz do backend já usado; nenhuma transferência GPU→CPU de campos foi adicionada. Foram alterados apenas os dois dataclasses/parsers, o encaminhamento pelos drivers, metadados e o ponto de taxa exclusivo da chuva. A instrumentação passiva anterior não foi alterada.

## Prova de neutralidade e testes

A função anterior foi congelada antes da edição em [`reference_rain_evaporation_legacy.py`](../tests/reference_rain_evaporation_legacy.py), com SHA256 do arquivo original no cabeçalho. Os auxiliares físicos que a referência importa permaneceram inalterados, inclusive constantes, termodinâmica, distribuição e ventilação.

[`test_rain_evaporation_factor.py`](../tests/test_rain_evaporation_factor.py) verifica:

- f=0: ausência de transferência de chuva;
- f=0.5 e f=1.1: proporcionalidade ao baseline sem limitadores ativos, com tolerância relativa de 3×10⁻¹⁶ para arredondamento da ordem de multiplicação;
- f=1: igualdade bit a bit da transferência com a função anterior;
- todas as demais funções de `PROCESS_ORDER`: transferências idênticas, no mesmo estado fixo, ao variar apenas f, incluindo casos com processos de gelo efetivamente ativos;
- três passos exclusivamente microfísicos com todos os processos: igualdade bit a bit de T e das sete espécies, comparando a rotina nova com a anterior, em CPU e GPU;
- preservação dos limites de água disponível e saturação, inclusive com fator grande de teste;
- validação, leitura de configurações antigas, serialização e encaminhamento do fator nos caminhos normal e de baixa memória, sem integrar tempestades nesses testes de configuração.

**62 testes passaram**, incluindo os 20 novos casos parametrizados, testes existentes de microfísica/acoplamento/termodinâmica e cinco testes da instrumentação. A suíte da instrumentação inclui passos curtos de verificação em malhas de teste; não são execuções da matriz causal nem uma tempestade completa. GPU/CuPy foi exercitada, sem skip nos testes novos.

Na avaliação local de **23.808 estados arquivados**, a transferência com f=1 foi idêntica à referência anterior: **erro máximo 0; RMS 0**. A neutralidade foi demonstrada para as rotinas e integrações microfísicas testadas. Não se afirma que uma nova tempestade completa de 3900 s foi executada para comparação.

## Sensibilidade local em estados reais da sequência

O script [`audit_rain_evaporation_control.py`](../scripts/audit_rain_evaporation_control.py) lê a sequência em modo somente leitura. Seleciona 31 saídas de aproximadamente 2400–3300 s, até 2 km, com 512 células subsaturadas contendo chuva e 256 células de cobertura geral por quadro. Os índices são determinísticos e uniformemente espaçados em cada estrato; pode haver sobreposição entre estratos. As estatísticas abaixo usam somente o estrato ativo: 15.872 entradas. Não são médias ponderadas da tempestade nem uma amostra aleatória independente.

Faixas amostradas: T=284,41–300,26 K; P=79,50–99,55 kPa; ρ=0,9525–1,1582 kg/m³; qv=0,00527–0,01512 kg/kg; qr=0–0,00728 kg/kg; z=39,89–1974,38 m. Os testes sintéticos de gelo cobrem adicionalmente 250–300 K; não se extrapolou a amostra quente para validar gelo.

Cada avaliação aplica **apenas uma transferência local de chuva, com Δt=0,5 s**, a partir do mesmo estado fixo. Esse passo fica dentro da faixa nativa arquivada de 0,4291–0,5904 s nessa janela. Não há transporte, movimento de parcelas, nova projeção, sedimentação ou integração de tempestade. Os estados são amostras sincronizadas pós-passo, não o estado intermediário exato que a rotina recebeu originalmente.

Calculam-se qv_new=qv+dq, qr_new=qr−dq, T_new=T−Lv/cp·dq, θ_new a partir de T_new e da mesma P, e θv_new=θ_new(1+0,61qv_new−Σq_cond,new). A flutuabilidade potencial usa os mesmos perfis de referência e gravidade do baseline: B=9,81[(θ−θ0)/θ0+0,61(qv−qv0)−Σq_cond]. “Potencial” significa a mudança que essa transferência produziria em B, antes de qualquer resposta dinâmica. Não é uma aceleração realizada por nova simulação.

| Mediana local em 0,5 s | f=0 | f=0,5 | f=0,9 | f=1 | f=1,1 |
|---|---:|---:|---:|---:|---:|
| Taxa efetiva (10⁻⁸ kg kg⁻¹ s⁻¹) | 0 | 2,2055 | 3,9699 | 4,4110 | 4,8521 |
| Δqv (10⁻⁸ kg/kg) | 0 | 1,1027 | 1,9849 | 2,2055 | 2,4260 |
| ΔT (10⁻⁵ K) | 0 | −2,7442 | −4,9396 | −5,4885 | −6,0373 |
| Δθ (10⁻⁵ K) | 0 | −2,8201 | −5,0762 | −5,6403 | −6,2043 |
| Δθv (10⁻⁵ K) | 0 | −2,3052 | −4,1494 | −4,6105 | −5,0715 |
| ΔB (10⁻⁶ m/s²) | 0 | −0,7434 | −1,3381 | −1,4868 | −1,6355 |
| Fração de entradas limitada | — | 0,743% | 1,720% | 1,896% | 2,092% |

Em f=1, o intervalo P05–P95 da taxa efetiva é 2,96×10⁻¹¹–1,39×10⁻⁶ kg kg⁻¹ s⁻¹. O intervalo P05–P95 de ΔT é −0,001734 a −3,68×10⁻⁸ K. Logo, a mediana pequena não representa os maiores efeitos locais.

Comparado a f=1, f=1,1 produz medianas das diferenças de aproximadamente −5,49×10⁻⁶ K em T e −1,49×10⁻⁷ m/s² em B nesse passo. f=0,9 produz diferenças opostas aproximadamente simétricas. A taxa bruta é escalada exatamente pela definição; a transferência limitada pode não ser proporcional, e θv contém um produto entre variáveis que também mudam.

**Esses resultados não escolhem ±10%.** Eles mostram uma resposta local pequena e quase linear em grande parte da amostra, com saturação dos limitadores em uma fração dos estados. Não estabelecem a magnitude integrada da perturbação do cold pool, nem sua resposta monotônica, nem efeitos na geometria do updraft. Não se pode multiplicar essas diferenças locais por centenas de segundos e tratar o resultado como previsão de uma tempestade.

## Influências diretas, indiretas e limites

Diretamente, só R_evap,rain recebe f. A transferência conserva água entre qr e qv e aplica o calor latente correspondente. A redução da carga de chuva e o aumento de vapor tendem a elevar B; o resfriamento tende a reduzi-la. O diagnóstico local inclui os três efeitos e encontra redução líquida de B na amostra ativa para f>1.

Indiretamente, os novos T, qv e qr podem modificar condensação, gelo, precipitação, sedimentação e circulação nos passos seguintes, apesar de suas equações permanecerem iguais. **Invariância dos outros processos no mesmo estado não significa invariância de sua evolução depois de uma intervenção.** A evaporação de nuvem, sublimação e derretimento continuam capazes de resfriar o ar quando f=0.

O novo fator não identifica exclusivamente “intensidade do cold pool”, não altera uma única componente de B e não fornece por si só uma conclusão causal sobre tilting ou manutenção de ζ. Ele identifica uma intervenção contínua em um processo de formação do cold pool, cuja mediação e efeitos concorrentes precisarão ser medidos em eventual etapa futura autorizada.

## Preservação e reprodução

O baseline, seus diagnósticos, configurações e metadados não foram regravados. O registro inicial de preservação contém o commit `097d516409285efad047012f80f8b2b8c615f71b`, o hash canônico da configuração, hashes dos fontes antes da edição e dos artefatos antigos. O baseline já continha alterações não commitadas de instrumentação; por isso o commit isolado não descreve toda a revisão. O manifesto final registra hashes dos fontes novos e dos arquivos entregues.

Resultados em [`outputs/rain_evaporation_control_audit`](../outputs/rain_evaporation_control_audit): `sampled_states.csv`, `local_responses.csv`, `local_statistics.csv`, `local_summary.json`, manifests e resultado de testes. [Figura da sensibilidade local](../outputs/rain_evaporation_control_audit/local_sensitivity.png).

```powershell
python scripts/audit_rain_evaporation_control.py outputs/diagnostic_sequence_20260905/sequence.h5 --out outputs/rain_evaporation_control_audit
python -m pytest tests/test_rain_evaporation_factor.py tests/test_microphysics.py tests/test_flow_microphysics_coupling.py tests/test_precip_microphysics_thermo_vectorization.py tests/test_diagnostic_capture.py -q
```

**Resposta final: SIM, COM LIMITAÇÕES.** Existe agora um controle contínuo e diretamente exclusivo da taxa de evaporação de chuva, numericamente neutro em f=1 nos testes realizados. Isso permite preparar um teste causal desse processo, mas não transforma evaporação de chuva em uma intervenção exclusiva sobre o cold pool. Nenhum caso causal foi executado e nenhum valor não padrão foi escolhido para uma futura tempestade.
