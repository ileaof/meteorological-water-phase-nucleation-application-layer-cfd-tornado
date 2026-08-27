# Prompt — Adicionar capacidade de simulação de tornadogênese

> Cole este prompt em um agente de código (Claude Code / Cowork) com o repositório
> `met_h2o_nucleation_cfd` aberto como contexto. Ele descreve **o quê**, **as
> restrições** e **os critérios de aceitação** — deixe o agente propor o plano
> detalhado e trabalhar de forma incremental, testando a cada etapa.

---

## Papel e objetivo

Você é um engenheiro sênior de CFD atmosférico. Trabalhando **dentro deste
repositório**, adicione a capacidade de **simular tornadogênese idealizada**: uma
supercélula que, a partir de um *sounding* com cisalhamento, desenvolve um
mesociclone de níveis médios e, em seguida, **rotação de baixo nível** perto da
superfície — o proxy físico de um tornado.

Seja honesto sobre o escopo: **isto é simulação idealizada, não previsão
operacional**. Nada de assimilação de dados, evento real ou condição inicial
observada nesta fase. O produto é um núcleo dinâmico *rotacional* validado contra
resultados clássicos de supercélula, reaproveitando a física já existente no repo.

## Restrições rígidas (não negociáveis)

1. **NÃO** modifique `src/met_water_nucleation/_engine/` — é o kernel de nucleação
   guardado por integridade (SHA-256). Após todo o trabalho, `--validate` (testes
   [1]–[21] + checksums) **deve continuar verde**.
2. **NÃO** altere o comportamento científico do pacote `meteorological_flow`
   atual nem quebre seus testes. Ele é o modelo "de demonstração" citável e deve
   permanecer reproduzível como está.
3. **Fork, não reescrita destrutiva.** Crie um **novo pacote** —
   `src/storm_dynamics/` (sugestão de nome; confirme comigo se preferir outro) —
   que *reutiliza* como biblioteca de física os módulos já existentes:
   `thermodynamics`, `precip_microphysics`, `base_state`, `grid`, `io`,
   partes de `diagnostics` e o kernel de nucleação. Reescreva apenas o **núcleo
   dinâmico**.
4. Mantenha o padrão do repo: `pyproject.toml` (adicionar o pacote em
   `packages.find`), testes em `tests/`, exemplo runnable em `examples/`, docs em
   `docs/`, config declarativa em `configs/`, e um `handoff.md` do novo pacote
   documentando decisões, estado e trabalho restante — no estilo do handoff atual.
5. Conservação preservada: água e energia devem continuar fechando aos mesmos
   níveis de erro dos budgets atuais; a continuidade de massa (resíduo da
   projeção) deve permanecer ~0.

## O que precisa ser construído (o núcleo do trabalho)

A razão de o modelo atual não poder girar está no `meteorological_flow._step`:
a advecção de momento está desligada, e há um *drag* de Rayleigh + teto de
velocidade que suprimiriam qualquer vórtice. O novo núcleo corrige isso. Em ordem
de habilitação:

### 1. Advecção de momento conservativa (o passo que habilita tudo)
Implemente o transporte advectivo completo do momento — `(u·∇)u` — em **flux-form
conservativo na grade staggered**, com reconstrução de 2ª ordem (MUSCL/minmod,
reaproveitando o estilo de `advection.advect_center_massflux`). Sem os termos de
**tilting** e **stretching** de vorticidade não há rotação; este item é
pré-requisito de todos os outros.

### 2. Força de Coriolis (plano-f)
Adicione Coriolis com `f` configurável (plano-f; latitude típica de *Tornado
Alley*, ~35–40°N). Documente que é plano-f idealizado.

### 3. Fechamento de turbulência LES (substituindo Rayleigh + clip)
Substitua o `gamma_damp` de Rayleigh e o `clip`/teto de 120 m/s por um
**fechamento subgrid LES** (Smagorinsky ou TKE-1.5 tipo Deardorff). A dissipação
passa a ser fisicamente motivada; nada de saturar velocidade artificialmente
(qualquer limite remanescente deve ser só um *guard* numérico extremo,
documentado).

### 4. Camada-limite / arrasto de superfície
No fundo do domínio, implemente uma **lei de arrasto bulk** (drag law) em vez de
apenas `free_slip`/`no_slip`. O atrito de superfície é essencial para a rotação de
baixo nível e a região de *corner flow*.

### 5. Cisalhamento direcional (hodógrafa curva)
O `base_state.weisman_klemp` atual só impõe cisalhamento unidirecional (`u_shear`).
Estenda para **hodógrafa curva** (quarter-circle / Klemp) populando `v0(z)` além de
`u0(z)`, gerando helicidade relativa à tempestade — condição necessária para
supercélula e rotação de baixo nível. A `BaseState` já tem os campos `u0`/`v0`.

### 6. Reuso da física existente
Acople a microfísica `precip_microphysics` já pronta (o resfriamento evaporativo /
cold pool é a fonte de vorticidade baroclínica da tornadogênese) e, opcionalmente,
o kernel de nucleação, exatamente como o `meteorological_flow` já faz. Não
reimplemente termodinâmica nem microfísica.

### 7. Diagnósticos de rotação
Adicione (estendendo `diagnostics`): vorticidade 3D, **vorticidade vertical ζ**,
*updraft helicity*, *storm-relative helicity* (SRH 0–1 km e 0–3 km), cisalhamento
0–1 km e 0–6 km, e rastreamento do **máximo de ζ perto da superfície** (o
indicador de tornadogênese). Esses são os números que provam que o modelo girou.

## Marcos (entregar e verificar um de cada vez)

- **M1 — Supercélula rotativa.** Sob cisalhamento unidirecional, o warm bubble deve
  produzir *storm splitting* (células defletidas à esquerda e à direita) e um
  mesociclone de níveis médios (ζ acima de um limiar em ~3–6 km). Este é o
  resultado clássico de Klemp–Wilhelmson / Weisman–Klemp e é o teste de sanidade
  de que a dinâmica de rotação funciona.
- **M2 — Rotação de baixo nível (proxy de tornadogênese).** Sob hodógrafa curva +
  arrasto de superfície + cold pool ativo, desenvolver ζ intenso perto da
  superfície na interface do *forward-flank downdraft*. Este é o alvo do projeto.
- **M3 — (stretch) Vórtice fino.** Refinamento aninhado / AMR para resolver o vórtice
  em escala de ~10–100 m. Provavelmente um projeto à parte; apenas deixe o núcleo
  preparado (não precisa entregar agora).

## Requisitos numéricos e de estabilidade

- CFL adaptativo como no núcleo atual, agora incluindo o momento advectado.
- Conservação de água/energia nos mesmos níveis atuais; resíduo de continuidade ~0
  vindo da projeção (não dos limitadores).
- Use o núcleo **anelástico** (já esboçado no repo) para a coluna profunda de
  tempestade — Boussinesq-esticado não é adequado a 10–12 km.
- Sem caps de velocidade não-físicos.

## Critérios de aceitação (definição de "pronto")

1. `--validate` verde e **todos os testes atuais passando** (o `_engine` e o
   `meteorological_flow` intactos).
2. Novos testes de regressão para M1 e M2: *storm splitting* reproduzido;
   mesociclone com ζ acima do limiar em níveis médios; rotação de baixo nível
   presente no cenário de hodógrafa curva; budgets de conservação dentro da
   tolerância.
3. Um exemplo runnable (`examples/supercell_tornadogenesis.py`) que roda em grade
   de demonstração e emite os diagnósticos de rotação (ζ, SRH, updraft helicity).
4. `docs/storm_dynamics_guide.md` explicando o modelo, o que ele pode e **o que
   NÃO pode** afirmar (idealizado, resolução, sem previsão de evento real), com as
   referências (Klemp & Wilhelmson 1978; Weisman & Klemp 1982; Rotunno & Klemp
   1985; núcleo não-hidrostático estilo CM1 de Bryan & Fritsch como referência).
5. `handoff.md` do novo pacote com decisões travadas, estado dos marcos e trabalho
   restante.

## Processo

Trabalhe **incrementalmente e com verificação a cada passo** — não escreva tudo de
uma vez. Antes de codar, apresente um plano curto: estrutura do novo pacote, quais
módulos existentes serão reusados como estão, e a ordem dos itens 1–7 acima. Rode
os testes após cada item. Ao terminar cada marco, mostre os diagnósticos de rotação
que provam o resultado.

## Não-objetivos (deixe explícito na doc)

Assimilação de dados, condições iniciais/contorno de eventos reais, previsão
operacional, e verificação contra observações **não** fazem parte desta fase. O
entregável é um simulador de **dinâmica de tempestade rotacional idealizada**
sobre a base de física já validada do repositório.
