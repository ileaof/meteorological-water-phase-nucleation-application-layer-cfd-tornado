# Validação do propagador discreto de proveniência de vorticidade

Data: 8 de setembro de 2026.

## Resposta

**O novo propagador de proveniência é discretamente consistente, estável e suficientemente bem condicionado para uma execução longa? SIM, dentro da definição explícita de partição adotada.**

O método aprovado é uma **partição linear do operando advectado com velocidade advectante e ramos MUSCL congelados a partir do estado total**. Ele recompõe exatamente o mapa advectivo nativo, salvo arredondamento, mantém condicionamento controlado por milhares de passos reduzidos e passou por 150 s da tempestade madura com prognósticos ON/OFF idênticos bit a bit.

Ele não é o Jacobiano completo do operador quadrático de Burgers. Essa diferença é matemática e necessária: para uma perturbação de escala, `J_M(U)U = 2M(U)`, enquanto uma decomposição cujas parcelas somam o próprio incremento exige `L_U(U)=M(U)`. Os rótulos devem ser interpretados como uma convenção aditiva de proveniência com advector comum, não como derivadas completas da solução em relação a cada fonte.

Não foram executados nesta etapa 0–3300 s, WEAK-EVAP, CONTROL ou STRONG-EVAP completos. Nenhum campo de proveniência foi usado para atribuição física.

## 1. Causa matemática da falha anterior

O traçador rejeitado avançava diretamente cada `ω_j` por Euler explícito e diferenças centradas:

```text
∂ω_j/∂t = −u·∇ω_j + (ω_j·∇)u − ω_j ∇·u.
```

Para advecção uniforme 1D, Euler centrado tem fator

```text
G = 1 − i C sin(kΔx),    |G| = sqrt(1 + C² sin²(kΔx)) > 1.
```

Esse propagador não reproduzia o mapa MUSCL em fluxo usado pelo momento. As parcelas cresciam secularmente e o `advection_remainder` adquiria quase a mesma norma e direção oposta. Em 2981 s, `κω≈3,98×10⁴` e `κT≈2,29×10⁵`. O fechamento algébrico escondia uma decomposição inutilizável.

No teste de referência atual, 300 passos do esquema antigo amplificaram a norma do modo de Fourier por `2,760642×10⁷`. O novo mapa permaneceu limitado.

## 2. Operador MUSCL real

O código auditado está em `src/storm_dynamics/momentum.py`. Para cada componente escalonada `q∈{u,v,w}`, o estágio explícito é

```text
qⁿ⁺¹ = qⁿ + Δt [−Dₓ(Aₓ Rₓq) − Dᵧ(Aᵧ Rᵧq) − D_z(A_z R_zq)].
```

`D` é a diferença conservativa dos fluxos no volume de controle da própria face. `A` é a velocidade transportadora interpolada para a face/canto do fluxo, e `R` é a reconstrução upwind MUSCL.

Para uma direção genérica, o slope nativo é

```text
d⁻ = q_i − q_{i−1}
d⁺ = q_{i+1} − q_i
s_i = 0                              se d⁻d⁺ ≤ 0
s_i = d⁻                             se |d⁻| < |d⁺|
s_i = d⁺                             caso contrário.
```

Na face:

```text
q_L = q_minus + 0,5 s_minus
q_R = q_plus  − 0,5 s_plus
q_face = q_L se A>0; q_R caso contrário
F = A q_face.
```

Não há clipping de velocidade nesse operador.

### Geometria C-grid

- Para `u`, o fluxo normal x usa `Uc=0,5(u_I+u_{I+1})` nos centros; os fluxos y e z interpolam `v` e `w` para os cantos do volume de controle de `u`.
- Para `v`, o fluxo normal y usa `Vc=0,5(v_J+v_{J+1})`; x e z são o espelho escalonado.
- Para `w`, o fluxo normal z usa `Wc=0,5(w_K+w_{K+1})`; `u` e `v` são interpolados para as faces z.
- Em z estendido, a divergência de `u` e `v` usa `dz_c`; a de `w` usa a distância entre centros `dzc_f`.
- Laterais periódicas usam `roll` e identificam as faces duplicadas. No modo não periódico, os fluxos externos são nulos e slopes nas bordas são zero.
- Z nunca é periódico. Os fluxos verticais nas paredes são zero e a tendência de `w` no solo é fixada em zero.

Assim, `Uⁿ⁺¹=M(Uⁿ)` é exatamente `U + Δt momentum_advection_tendency(U)`, com todos os fluxos calculados no posicionamento nativo de cada componente.

## 3. Formulação escolhida

Cada fonte é armazenada primeiro como velocidade escalonada:

```text
U = Σ_j U_j.
```

Em cada estágio advectivo, o estado total define:

- velocidades transportadoras;
- sinal upwind;
- ramo minmod;
- métricas e tratamento de contorno.

Essas decisões são congeladas e o mesmo mapa linear `L_U` é aplicado a todos os `U_j`:

```text
U_jⁿ⁺¹ = U_jⁿ + Δt L_U(U_jⁿ)
Σ_j L_U(U_j) = L_U(Σ_j U_j) = M_adv(U).
```

Depois de cada estágio, diagnostica-se `ω_j=curl(U_j)` usando as mesmas posições C-grid. Isso é superior ao transporte direto de `ω_j`: evita misturar centros e faces e herda exatamente a geometria do fluxo prognóstico.

## 4. Tratamento rigoroso do limiter

O ramo é sempre o escolhido pelo estado total naquele estágio:

- mudança de sinal ou qualquer slope zero: ramo zero, pois `d⁻d⁺≤0`;
- empate de módulos não nulos: ramo forward, reproduzindo o `else` nativo;
- velocidade transportadora positiva: estado esquerdo;
- velocidade zero: o código seleciona o estado direito, mas o fluxo é exatamente zero;
- borda rígida: slope zero;
- não diferenciabilidade: usa-se a derivada direcional do ramo efetivamente selecionado, sem tentar atravessar a descontinuidade do limiter.

O estudo de ε, longe das superfícies de troca, comparou `_frozen_slope` com diferenças finitas do `minmod` nativo. O erro máximo foi `6,20×10⁻¹³` em `ε=10⁻²`; abaixo de aproximadamente `10⁻⁶`, o erro cresce como esperado por cancelamento de ponto flutuante, chegando a `3,88×10⁻⁷` em `10⁻⁸` e `5,70×10⁻⁶` em `10⁻⁹`. Não houve troca de ramo escondida no residual.

## 5. Jacobiano completo versus partição aditiva

Foi executada a diferença finita solicitada sobre um campo não uniforme:

```text
[M(U+εδU)−M(U)]/ε.
```

Para `δU=U`, a homogeneidade quadrática fornece `J_M(U)U=2M(U)`. O erro relativo para `2M` caiu de `1,00×10⁻²` em `ε=10⁻²` para `1,00×10⁻⁶` em `ε=10⁻⁶`; em `10⁻⁸`, arredondamento elevou-o a `2,42×10⁻⁶`. A diferença entre o Jacobiano completo e a partição congelada convergiu para 50%.

Portanto, usar `J_M` integralmente em cada origem duplicaria a contribuição bilinear e impediria `ΣU_j=U`. O método implementado é a decomposição equivalente `F(U,ΣU_j)=ΣF(U,U_j)`. Ela é exata, estável e reprodutível, mas não é única: outra convenção poderia repartir o papel da velocidade transportadora entre rótulos.

## 6. Testes 1D

Foram propagados onda suave, pulso, gradiente, mudança de sinal e campo com limiter ativo por 1.600 passos, equivalentes a 8,75 travessias periódicas.

| caso | máximo erro de soma | razão amplitude final/inicial | `κ` máximo |
|---|---:|---:|---:|
| onda suave | 1,20×10⁻¹⁴ | 0,9901 | 1,00384 |
| pulso | 8,10×10⁻¹⁵ | 0,9877 | 1,00416 |
| gradiente periódico | 5,16×10⁻¹⁵ | 0,8325 | 1,00631 |
| mudança de sinal | 1,29×10⁻¹⁴ | 1,0000 | 1,00416 |
| limiter ativo | 2,49×10⁻¹⁴ | 0,7628 | 1,00782 |

Não houve crescimento explosivo. A redução de amplitude nos casos descontínuos/alta frequência é a dissipação do MUSCL limitado, comum ao total e às parcelas.

## 7. Testes multidimensionais

Translação 3D, shear, rotação e deformação periódica incompressível foram executados por 250 passos. Todos permaneceram finitos, sem explosão de norma, e recompuseram o total com erro relativo inferior a `2×10⁻¹²`.

O primeiro campo de deformação usado durante o desenvolvimento era descontínuo na fronteira periódica e tinha compressão efetiva; seu crescimento não foi atribuído ao propagador. O teste definitivo usa strain periódico analiticamente não divergente.

## 8. Teste C-grid

Um campo escalonado foi decomposto em três parcelas não proporcionais. Foram verificados simultaneamente:

```text
momentum_advection_tendency(U) = Σ_j L_U(U_j)
curl(Σ_j U_j) = Σ_j curl(U_j).
```

O primeiro é bit a bit idêntico quando o único rótulo é o total e fecha até arredondamento para a soma não trivial. O segundo fecha com tolerâncias relativas/absolutas de `2×10⁻¹³`, usando o mesmo operador de curl da saída diagnóstica.

## 9. Fontes isoladas

Testes separados para buoyancy, LES, surface drag, Coriolis e projection injetaram um incremento nativo conhecido. Em todos os casos:

1. somente o rótulo correspondente recebeu o incremento;
2. os demais rótulos permaneceram exatamente zero;
3. o campo foi alterado pelo estágio advectivo posterior;
4. o fechamento relativo de `ω` permaneceu abaixo de `10⁻¹¹`.

Na execução real, boundary e other seguem o mesmo mecanismo de captura dos incrementos efetivamente aplicados. O `advection_remainder` recebe somente a diferença mensurada entre o incremento nativo e a soma congelada; nenhum coeficiente é ajustado para fazê-lo pequeno.

## 10. Estabilidade de longa memória reduzida

O teste C-grid reduzido executou 2.500 passos:

| métrica | valor |
|---|---:|
| `κω` inicial | 1,000000059 |
| `κω` final | 1,000000040 |
| `κω` máximo | 1,000000059 |
| máximo erro absoluto de curl | 7,99×10⁻¹⁸ s⁻¹ |

Não há tendência exponencial ou monotônica de crescimento. A comparação temporal CFL 0,10 versus 0,05, no mesmo tempo final, produziu diferenças RMS relativas de 1,783% no total, 1,788% no rótulo 1 e 1,797% no rótulo 2. A sensibilidade adicional das parcelas em relação à própria solução foi inferior a 0,02 ponto percentual.

## 11. Tempestade madura, sem execução completa

O snapshot arquivado de 2790,253490 s foi avançado por 335 passos nativos até 2940,145878 s na malha `120×120×48`, em GPU. Foram executadas lado a lado uma cópia sem instrumento e outra instrumentada.

- `u`, `v`, `w`, `p`, `theta`, `qv`, `ql`, `qi`, `qr`, `qs`, `qg` e `qh`: idênticos bit a bit em todos os passos verificados;
- maior diferença prognóstica absoluta: zero;
- `RMS(Rω)/RMS(ω)=3,31×10⁻¹⁵`;
- `max|Rω|=4,23×10⁻¹⁶ s⁻¹` no final e `4,82×10⁻¹⁶ s⁻¹` no máximo da janela;
- `RMS(RT)/RMS(T)=3,39×10⁻¹⁵`;
- `max|RT|=1,21×10⁻¹⁸ s⁻²` no final e `1,33×10⁻¹⁸ s⁻²` no máximo da janela;
- `κω=1,2673` e `κT=1,6551` no final;
- `||Radv,ω||/||ω||=3,80×10⁻¹⁵`;
- `||Radv,ω||/Σ||ω_j||=3,00×10⁻¹⁵`;
- `||Radv,T||/||T||=3,77×10⁻¹⁵`;
- `||Radv,T||/Σ||T_j||=2,28×10⁻¹⁵`.

Todos os limites de fechamento originais foram satisfeitos por pelo menos três ordens de grandeza; na prática, por quatro a cinco.

## 12. Condicionamento e critérios empíricos

Não foi escolhido um teto de `κ` para aprovar retrospectivamente. Os valores de referência observados foram:

- partições sintéticas bem condicionadas: `κ≈1,000–1,008`;
- tempestade madura com fontes físicas distintas: `κω=1,267`, `κT=1,655` após 150 s;
- instrumento antigo inválido: `κω≈3,98×10⁴`, `κT≈2,29×10⁵`.

O método novo não mostra crescimento secular nos testes reduzidos e permanece várias ordens de grandeza abaixo do regime de cancelamento anterior no segmento maduro. Uma execução longa futura deve continuar monitorando `κ(t)` sem transformar estes valores observados em limite universal.

## 13. Neutralidade e custo

No teste maduro:

- tempo sem instrumento: 171,77 s;
- tempo instrumentado: 607,38 s;
- razão: 3,536;
- memória persistente dos nove rótulos: 184.550.400 bytes, aproximadamente 176 MiB;
- pico medido do dispositivo com duas simulações coexistentes: 2.359.820.288 bytes;
- arquivo HDF5 da janela: 215.807.158 bytes.

O custo de produção com apenas uma simulação será menor que o pico pareado, mas o overhead temporal continua material.

## 14. Limitações interpretativas

- A partição frozen-advector é uma convenção aditiva coerente, não uma decomposição única do operador bilinear.
- O Jacobiano completo é útil para sensibilidade, mas não satisfaz a soma exigida para estes rótulos.
- Empates e mudanças de ramo do limiter não têm derivada clássica; a regra frozen-branch é explícita e determinística.
- O teste de 150 s não substitui o monitoramento contínuo em uma futura execução de memória completa.
- Aprovação matemática do instrumento não autoriza interpretar arquivos produzidos por versões anteriores nem autoriza uma nova matriz causal.

## 15. Artefatos e reprodução

- Implementação: `src/storm_dynamics/vorticity_provenance.py`;
- testes básicos: `tests/test_vorticity_provenance.py`;
- testes estendidos: `tests/test_vorticity_provenance_propagator.py`;
- validação reduzida: `scripts/validate_vorticity_provenance_propagator.py`;
- resultados reduzidos: `outputs/vorticity_provenance_propagator_validation`;
- segmento maduro: `outputs/vorticity_provenance_validation_mature_150s/summary.json`.

Comandos:

```powershell
python -m pytest tests/test_vorticity_provenance.py tests/test_vorticity_provenance_propagator.py -q
python scripts/validate_vorticity_provenance_propagator.py
python scripts/validate_vorticity_provenance_mature.py --snapshot 93 --steps 335 --device gpu --out outputs/vorticity_provenance_validation_mature_150s
```

## Decisão final

**APROVADO PARA EXECUÇÃO LONGA**

A aprovação vale para o instrumento v4 e exige manter o monitoramento contínuo de fechamento, ambos os índices de condicionamento e as duas normalizações de `advection_remainder`. Ela não autoriza por si só iniciar uma execução longa ou interpretar experimentos anteriores.
