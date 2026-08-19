# Methods for Ranking Nearby Stars as Solar Gravitational Lens SETI Targets

## Purpose

This note describes heuristic methods for ranking nearby stars whose **solar gravitational-lens (SGL) focal lines** could host equipment associated with a hypothetical interstellar communications network.

The motivating network consists of autonomous probes that expand from star to star, establish bidirectional SGL links, and add redundant connections. Before running full network-growth simulations, an archival search can examine dozens of candidate focal lines relatively cheaply. The near-term objective is therefore not to identify one uniquely “correct” target list. It is to build a **high-recall portfolio** that covers several plausible network architectures.

The central recommendation is:

> Keep the nearest-star list as a protected baseline, then add stars favored by engineering quality, natural geometric adjacency, hub or gateway value, historical geometry, and archive searchability.

A robust target list should be the union of several rankings rather than the top entries from one heavily assumption-dependent score.

---

## 1. What is being ranked?

A single catalog row called a “star” can represent several different decisions. These should be separated.

### 1.1 Stellar-system edge probability

How likely is a hypothetical network to contain an active link between the Solar System and a neighboring **stellar system**?

This is primarily a graph and expansion question. Relevant factors include distance, geometric adjacency, expansion history, redundancy value, and whether the system is an attractive relay or replication site.

### 1.2 Lens-component probability

If the neighboring system contains multiple stars, which component would be used as the gravitational lens?

For a component \(c\) in system \(S\), a useful conceptual decomposition is

\[
P(\odot \leftrightarrow c)
=
P(\odot \leftrightarrow S)
P(c\text{ is selected as the lens}\mid \odot \leftrightarrow S).
\]

Close binary components should generally count as one industrial or graph node, but they may generate separate component-specific focal-line tracks for an archival search.

### 1.3 Archive-search priority

Even a physically plausible link may be a poor archival target because its focal-line locus is badly covered, confused, poorly localized, or observable only under a narrow transmission geometry.

A useful decomposition is

\[
\begin{aligned}
Y_j \propto {} & P(E_j) \\
&\times P(E_j\text{ is active at the archive epoch}\mid E_j) \\
&\times P(\text{detectable signature}\mid E_j) \\
&\times Q_{\mathrm{archive},j},
\end{aligned}
\]

where \(E_j\) is a Solar-System–candidate-system edge.

The **network prior** and **archive yield** should be calculated separately and combined only at the end.

---

## 2. SGL geometry that affects the ranking

The Sun’s electromagnetic gravitational focal region begins at roughly 548 AU and continues outward as a focal half-line rather than terminating at one focal point.[^1]

For a lens star of mass \(M_\star\) and radius \(R_\star\), the approximate minimum focal distance is

\[
z_{\min,\star}
\simeq
\frac{R_\star^2 c^2}{4GM_\star}.
\]

This leads to a useful compactness proxy:

\[
\mathcal{C}_{\mathrm{focal}}
=
\frac{M_\star}{R_\star^2},
\qquad
z_{\min,\star}\propto \mathcal{C}_{\mathrm{focal}}^{-1}.
\]

Two qualifications are important:

1. **The target star does not change the Sun-side minimum focal distance.** Every candidate Sun–star link uses the Sun as the lens at our end. The candidate’s mass and radius affect the remote endpoint and therefore the economic plausibility of the whole bidirectional edge.
2. **A small minimum distance is optional, not mandatory.** A civilization can place equipment farther out. Large \(z_{\min}\) is restrictive; small \(z_{\min}\) merely expands the operating choices.

A better engineering score therefore evaluates an assumed operating distance \(z_{\mathrm{op}}\), not only \(z_{\min}\).

---

## 3. General ranking philosophy

### 3.1 Optimize for recall before precision

When dozens of focal lines can be checked, the cost of a false positive target is low relative to the cost of omitting a genuinely likely link. The first target set should therefore:

- protect the nearest systems from over-aggressive downranking;
- include candidates favored under several different civilization objectives;
- reserve a small wildcard allocation for unusual architectures;
- record why each target was included.

### 3.2 Avoid one universal weighted score

A single score can conceal structural assumptions. A civilization minimizing construction energy, one minimizing latency, and one maintaining an ancient historical graph can select substantially different neighbors.

Use several scenario-specific rankings, then aggregate them. Candidates recurring across incompatible scenarios are the strongest general-purpose targets. Candidates appearing only in one unusual scenario remain useful as explicitly labeled wildcards.

### 3.3 Treat uncertain observations honestly

A missing companion, planet, rotation period, or radial velocity is not the same as a favorable measurement. Distinguish:

- measured favorable properties;
- measured unfavorable properties;
- unknown properties;
- upper or lower limits;
- data-quality or multiplicity warnings.

---

## 4. Proximity and kinematic criteria

## 4.1 Current distance

Current Sun–star distance should remain the baseline and will probably dominate many rankings. Shorter edges generally imply:

- shorter probe travel;
- lower launch energy or faster activation;
- lower communication latency;
- lower interstellar propagation loss;
- greater likelihood of selection by local expansion rules.

Possible distance terms include

\[
S_d=d^{-\alpha}
\]

or a soft maximum-link penalty. Different \(\alpha\) values represent different assumptions about how strongly the civilization values short links.

Do not allow engineering penalties to remove all of the closest stars from the initial search portfolio. An advanced civilization may tolerate engineering inconvenience, and a persistent construction graph may retain an originally cheap parent edge even after better alternatives become available.

## 4.2 Construction-era and closest-approach distance

For an old or persistent network, present geometry may not be the relevant geometry. Calculate, over an assumed age window \(H\):

\[
d_{\min}(H)=\min_{-H<t<0}d_{\odot j}(t),
\]

and optionally

\[
f_k(H)=
\text{fraction of the interval during which }j
\text{ was among the Sun's }k\text{ nearest systems}.
\]

Useful parallel rankings are:

- **present-geometry:** favors a young or frequently rewired network;
- **construction-epoch:** samples a plausible time when the Sun was connected;
- **closest-approach:** favors opportunistic or historical links;
- **time-averaged:** favors long-lived economical edges.

## 4.3 Full relative motion

Use full three-dimensional relative velocity where possible. Proper motion alone can be misleading because it depends on distance and projection.

Useful quantities include:

- present relative position and velocity;
- change in Sun–star separation over the assumed service interval;
- time spent within an economic link radius;
- uncertainty from missing or imprecise radial velocity;
- persistence of natural-neighbor status over time.

## 4.4 Proper motion versus stationkeeping

A constant proper motion mainly sets a planned transverse velocity for a focal-line node. At Sun-centered distance \(z\), the line’s transverse speed scales approximately as

\[
v_{\perp}\simeq z\mu.
\]

The more relevant ongoing burden is angular acceleration and unpredictable or periodic motion:

\[
a_{\perp}\sim z\left\lVert\ddot{\hat{\mathbf n}}\right\rVert.
\]

Potential contributors include:

- perspective acceleration from radial motion;
- orbital motion in a multiple system;
- planetary reflex motion;
- unresolved companions;
- long-term gravitational accelerations;
- uncertainty in the ephemeris.

Accordingly, high proper motion should not be an automatic veto. It is more useful as a search-localization and edge-persistence parameter than as a direct fuel proxy.

## 4.5 Kinematic persistence score

Propagate the local catalog through a set of plausible construction epochs and recompute simple geometric graphs. Define

\[
K_j=
P\left[(\odot,j)\text{ is a natural-neighbor edge at the relevant epoch}\right].
\]

This is much cheaper than a network-growth simulation and can reveal stars that remain geometrically important over long periods.

---

## 5. Multiplicity and dynamical cleanliness

## 5.1 Do not use a binary yes/no penalty

The relevant engineering variable is the acceleration of the selected lens star caused by its companions. A rough proxy is

\[
a_{\mathrm{reflex}}
\sim
\sum_k \frac{GM_k}{a_k^2},
\]

with orbital phase, eccentricity, inclination, and projected geometry included when known.

A practical qualitative classification is:

| Architecture | Default treatment |
|---|---|
| Contact or very close binary | Strong negative |
| Close, comparable-mass binary | Strong negative |
| Moderate-period stellar binary | Moderate to strong negative |
| Hierarchical multiple with a distant outer companion | Weak to moderate negative |
| Very wide companion | Often weak; evaluate actual acceleration |
| Apparently single star | Positive, but retain unseen-companion uncertainty |

Close stellar companions can force the lens star to execute significant periodic reflex motion and can complicate the lens potential. Rapid rotation can also make the lens appreciably nonspherical.[^2]

## 5.2 Suggested dynamical-hostility metric

A normalized score could use

\[
Q_{\mathrm{dyn}}
=
-\log\left(1+\frac{a_{\mathrm{reflex,rms}}}{a_0}\right),
\]

where \(a_0\) represents a scenario-dependent tolerable acceleration.

When orbit information is incomplete, use broad categories and uncertainty ranges rather than false precision.

## 5.3 Planets and acceleration trends

Known planet count is a poor direct ranking variable because discovery completeness varies strongly by target. Prefer:

- measured astrometric acceleration;
- long-term radial-velocity trends;
- known massive-companion parameters;
- upper limits on unseen companions;
- non-single-star or excess-noise indicators.

Planetary systems may be positive for resources while negative for reflex dynamics. Keep those terms separate.

---

## 6. Lens-host engineering criteria

## 6.1 Minimum focal distance

Use \(R_\star^2/M_\star\) as a transparent first-order term. Strongly penalize very large values, especially for evolved giants. Give only a limited bonus to very small values because the civilization may choose to operate farther out.

A practical asymmetric scoring rule is:

- normal main-sequence values: modest differentiation;
- unusually compact main-sequence stars: small bonus;
- giants and supergiants: strong penalty;
- compact remnants: separate architecture-dependent category.

## 6.2 Operating-distance scenarios

Evaluate several possible operating rules:

| Policy | Example definition | Ranking consequence |
|---|---|---|
| Minimum-distance | \(z_{\mathrm{op}}=kz_{\min}\) | Strongly favors compact stars |
| Standard relay distance | \(z_{\mathrm{op}}=\max(z_{\min},z_{\mathrm{standard}})\) | Reduces differences among ordinary dwarfs |
| Low-gravity operation | Require \(GM/z_{\mathrm{op}}^2<a_{\max}\) | Can penalize massive compact objects |
| Starlight-powered | Require \(L/z_{\mathrm{op}}^2>F_{\min}\) | Favors luminous stars |
| Multi-objective | Minimize travel, gravity, power, background, and maintenance | Most flexible, most assumption-dependent |

This scenario separation is preferable to declaring one stellar class universally optimal.

## 6.3 Rotation and oblateness

Possible metrics include:

- rotation period;
- projected rotation speed \(v\sin i\);
- estimated fraction of breakup speed;
- estimated quadrupole or oblateness;
- evidence for time-variable surface asymmetry.

Rapid rotation is a plausible complexity and stationkeeping penalty, but not necessarily a hard exclusion for a very advanced network.

## 6.4 Stellar activity, corona, and wind

Activity can matter through plasma, radiation, variability, and natural background. Possible inputs include:

- X-ray to bolometric luminosity ratio;
- flare frequency and amplitude;
- chromospheric activity indicators;
- wind or mass-loss proxies;
- stellar age and rotation.

The weight should depend on the assumed communication wavelength. Coronal plasma is more troublesome for radio than for optical or near-infrared links, while stellar brightness and suppression requirements remain relevant at shorter wavelengths.[^1]

## 6.5 Power and logistics

Do not assume that the focal node must be powered by local starlight. Evaluate separate architectures:

- direct stellar collection at the node;
- beamed power from inner-system infrastructure;
- fission or fusion fuel transported outward;
- large long-life power stores;
- distributed collectors and a comparatively lightweight focal antenna.

For a direct-starlight architecture, a simple term is

\[
Q_{\mathrm{power}}
=
\log\left(\frac{L_\star}{z_{\mathrm{op}}^2}\right).
\]

For fuel- or beam-powered networks, the star’s luminosity may matter much less than the transport distance, total system resources, and reliability of the inner-system industrial base.

---

## 7. Longevity and stellar evolution

The useful quantity is not spectral class by itself but the probability of disruptive evolution during the intended service horizon:

\[
Q_{\mathrm{life}}
=
\min\left(1,
\frac{t_{\mathrm{stable,remaining}}}{T_{\mathrm{required}}}
\right).
\]

Once a star’s stable lifetime greatly exceeds the network’s service horizon, additional longevity should provide little extra score.

A broad default treatment is:

| Stellar category | Heuristic treatment |
|---|---|
| Old, isolated K dwarf | Strong conventional candidate |
| Old, inactive early or mid M dwarf | Strong candidate; activity and power caveats |
| Isolated G dwarf | Favorable or neutral |
| F or A dwarf | Soft to moderate penalty for lifetime, size, and often rotation |
| Subgiant or giant | Strong penalty for size and evolutionary change |
| Very young star | Penalty for activity and environmental instability |
| Isolated white dwarf | High-variance wildcard |
| Accreting white dwarf or compact interacting binary | Strong negative |

For most ordinary nearby main-sequence dwarfs, explosion risk is not the primary discriminator over modest network lifetimes. Longevity becomes most useful as an extreme-object filter.

## 7.1 White dwarfs as a separate hypothesis

Isolated white dwarfs have very small radii and substantial masses, producing extremely favorable formal focal minima. They are also dim and long-lived as remnants. Possible disadvantages include strong gravity at very small operating distances, low natural power, prior disruption of the planetary system, and limited local resources.

Do not mix them uncritically into the conventional ranking. Maintain a small **compact-lens architecture** basket instead.

---

## 8. Static geometric graph criteria

These methods estimate graph value without simulating autonomous expansion.

## 8.1 Mutual nearest-neighbor rank

For candidate \(j\), let:

- \(r_\odot(j)\) be its rank among the Sun’s nearest systems;
- \(r_j(\odot)\) be the Sun’s rank among the candidate’s nearest systems.

A simple mutuality score is

\[
S_{\mathrm{mutual}}(j)
=
\frac{1}{r_\odot(j)}
+
\frac{1}{r_j(\odot)}.
\]

A system moderately far from the Sun can still be a natural direct neighbor if the Sun is unusually close from its perspective.

## 8.2 Delaunay and Voronoi adjacency

In a three-dimensional Delaunay graph, two stellar systems are adjacent when their Voronoi cells share a face. This is an especially useful natural-edge heuristic because it identifies systems occupying distinct spatial sectors around the Sun, including some not at the very top of the raw distance list.

Possible refinements include:

- shared Voronoi-face area;
- stability of the adjacency under astrometric uncertainties;
- persistence under stellar motion;
- whether the edge also appears in stricter proximity graphs.

## 8.3 Gabriel and relative-neighborhood graphs

For a Sun–star candidate edge:

- **Gabriel condition:** no other eligible system lies inside the sphere whose diameter is the Sun–star segment.
- **Relative-neighborhood condition:** there is no third eligible system closer to both endpoints than the endpoints are to one another.

Edges surviving these tests are plausible in sparse, energy-conscious networks.

## 8.4 Natural-adjacency voting

Construct several static graph families and score recurrence:

- mutual \(k\)-nearest-neighbor graphs for several \(k\);
- maximum-distance graphs for several link horizons;
- Delaunay adjacency;
- Gabriel and relative-neighborhood graphs;
- Euclidean minimum-spanning trees;
- spanning trees augmented until bridges are removed;
- host-weighted versions of the above.

Define

\[
A_j
=
\sum_g w_g\,
\mathbf{1}[(\odot,j)\in g].
\]

A candidate repeatedly connected to the Sun across different sparse graph families is a robust geometric target.

## 8.5 Filtered adjacency for a picky network

To model a network that skips unfavorable systems:

1. Calculate a host-quality score.
2. Remove or strongly penalize the worst fraction of systems.
3. Recompute the proximity graphs.
4. Record which farther systems become direct Sun neighbors.

A selective-adjacency score can be written as

\[
A_j^{\mathrm{selective}}
=
\sum_q\sum_g
w_{qg}\,
\mathbf{1}[(\odot,j)\in g\text{ after eligibility cut }q].
\]

This is one of the best pre-simulation methods for finding plausible farther links that a nearest-20 list could miss.

---

## 9. Hub and gateway criteria

A star can be valuable in two different ways.

## 9.1 Dense hub

A dense hub has many cheap onward connections. Possible metrics include

\[
\rho_k(j)\propto \frac{k}{d_{j,k}^3},
\]

where \(d_{j,k}\) is the distance to the candidate’s \(k\)-th eligible neighboring system, and simple neighbor counts

\[
N_j(r)=\#\{u:d(j,u)<r\}.
\]

Calculate density using **stellar systems**, not individual components in a close multiple.

## 9.2 Gateway hub

A gateway provides access to a region poorly reached directly from the Sun. A simple incremental-reach score is

\[
I_j(r)
=
\#\{u:d(j,u)<r\ \land\ d(\odot,u)\ge r\}.
\]

Other gateway proxies include:

- weighted betweenness in a sparse geometric graph;
- stars reachable cheaply from the candidate but not from likely nearer Sun neighbors;
- access to a different spatial or angular sector;
- reduction in bridge vulnerability when the edge is added;
- large latency improvement relative to indirect routes.

## 9.3 Shortcut value

Construct a plausible local graph without the direct Sun–candidate edge. Let \(D_{\mathrm{alt}}\) be the shortest existing path length between them. Then

\[
S_{\mathrm{shortcut}}
=
\frac{D_{\mathrm{alt}}}{d_{\odot j}}.
\]

Large values identify direct edges that create major latency shortcuts. These candidates are particularly relevant to backbone-oriented networks.

## 9.4 Do not promote generic centrality blindly

A globally central nearby star is not automatically a likely direct Sun neighbor. Require at least one direct-edge justification:

- geometric adjacency;
- acceptable link distance;
- mutual-neighbor status;
- major shortcut value;
- large arrival-direction share;
- strong role after unfavorable intermediates are removed.

---

## 10. Construction-history proxies without a full simulation

## 10.1 Incoming-wave parent probability

If the mature network retains parent links, the Sun’s parent may reflect the direction from which the expansion wave arrived.

For each hypothetical arrival direction:

1. classify candidate systems on the upstream side of the Sun;
2. choose a preferred upstream parent using distance and host cost;
3. repeat across directions;
4. measure each candidate’s share of the sky.

Define

\[
W_j
=
\frac{\text{solid angle of arrival directions for which }j
\text{ is the preferred parent}}{4\pi}.
\]

This can be repeated with:

- all systems eligible;
- poor hosts removed;
- different maximum probe ranges;
- different distance and host-quality weights;
- present and historical stellar positions.

A candidate with a large \(W_j\) is a useful archival target even if it is not extremely close.

## 10.2 Persistent versus rewired networks

Maintain separate heuristic rankings for:

- **persistent construction graphs**, where old parent edges remain active;
- **mature rewired graphs**, where edges are replaced to fit current geometry and engineering cost.

Persistent rankings should emphasize construction-era distance, arrival direction, and historical proximity. Rewired rankings should emphasize current geometry, current host quality, and present static graph value.

---

## 11. Rank likely neighbor pairs as well as individual stars

The Sun is assumed to have at least two links. The value of one neighbor may depend on the other.

## 11.1 Cheap-cycle pairs

A pair \((j,k)\) may form an inexpensive local triangle. Favor:

- short Sun–\(j\) and Sun–\(k\) edges;
- short \(j\)–\(k\) separation;
- mutual natural-neighbor status;
- all three systems being dynamically acceptable.

A conceptual score is

\[
S_{jk}^{\mathrm{triangle}}
=
S_j+S_k
-\lambda_1(d_{\odot j}+d_{\odot k})
-\lambda_2d_{jk}.
\]

## 11.2 Independent-branch pairs

A resilience-oriented network may instead prefer neighbors serving different regions. Favor:

- low overlap between their onward neighborhoods;
- access to different geometric clusters or angular sectors;
- independent routes to the wider graph;
- low shared failure risk;
- one cheap local edge plus one farther gateway edge.

Stars recurring in strong pairs should be promoted even when their individual rank is only moderate.

---

## 12. Resource and industrial criteria

These are plausible but generally weaker and more observationally biased than geometry and dynamics.

Potential indicators include:

- known planets or debris reservoirs;
- total system mass;
- metallicity as a weak proxy for solid material;
- system age and dynamical stability;
- availability of multiple resource-bearing bodies;
- feasibility of beaming or transporting power to the focal region.

Cautions:

- planet and debris-disk catalogs are incomplete;
- a detected debris disk may indicate both raw material and collision or background problems;
- metallicity is not a direct measurement of accessible industrial inventory;
- known giant planets may help resources while worsening reflex motion.

These should receive low default weights or be placed in explicit resource-rich scenarios.

---

## 13. Archive-specific criteria

## 13.1 Focal-line localization

For each candidate, propagate the target direction and map the antipodal Sun-side locus over:

- the archive’s observation dates;
- an assumed range of heliocentric node distances;
- Earth’s orbital parallax;
- the star’s proper motion and radial motion;
- target astrometric covariance;
- component orbital motion in multiple systems.

The expected source is a track or search tube, not a static point. The first published radio search aimed at an Alpha Centauri SGL relay explicitly corrected for Earth-orbit parallax while searching beyond the Sun’s innermost focal distance.[^3]

## 13.2 Survey coverage

For each signature class, record:

- sky coverage;
- limiting flux or brightness;
- angular and spectral resolution;
- number of epochs;
- cadence and time baseline;
- masking and saturation;
- moving-object or transient completeness;
- whether the focal-track uncertainty is fully covered.

## 13.3 Confusion and natural backgrounds

Potential penalties include:

- stellar crowding;
- Galactic-plane background;
- diffuse infrared cirrus;
- nearby saturated stars;
- strong radio-frequency interference;
- natural variable or narrowband emitters;
- poorly characterized detector artifacts.

These alter search efficiency, not the physical probability of a link.

## 13.4 Signature-specific geometry

Maintain separate archive rankings for:

- radio leakage or local relay transmissions;
- tightly directed primary communication beams;
- optical or near-infrared laser leakage;
- thermal infrared emission;
- reflected light from an artifact;
- propulsion or stationkeeping signatures;
- anomalous moving objects.

Ecliptic or Earth-crossing alignment can be important for intercepting a narrow primary beam but much less important for thermal emission or direct imaging of the equipment.

---

## 14. Scenario-specific scoring profiles

For each candidate \(j\), calculate normalized subscores such as

\[
\mathbf{s}_j=
(D_j, A_j, A_j^{\mathrm{selective}}, H_j,
Q_{\mathrm{host},j}, K_j, W_j, Q_{\mathrm{archive},j}).
\]

Where:

- \(D\): proximity and link cost;
- \(A\): natural geometric adjacency;
- \(A^{\mathrm{selective}}\): adjacency after poor hosts are filtered;
- \(H\): hub or gateway value;
- \(Q_{\mathrm{host}}\): engineering quality;
- \(K\): kinematic persistence;
- \(W\): incoming-wave parent share;
- \(Q_{\mathrm{archive}}\): archive searchability.

Suggested qualitative weight profiles are:

| Network hypothesis | Distance | Host quality | Natural adjacency | Hub value | Historical geometry |
|---|---:|---:|---:|---:|---:|
| Inclusive persistent network | Very high | Low | High | Medium | High |
| Picky economical network | High | Very high | Very high after filtering | High | Medium |
| Latency backbone | Medium | High | Medium | Very high | Medium |
| Energy-frugal sparse network | Very high | High | Very high | Low to medium | Medium |
| Ancient mature rewired network | High at observation epoch | High | Very high | High | Very high |
| Compact-lens architecture | Medium | Dominated by compactness and dynamics | Medium | High | Medium |

Possible combination methods include:

- scenario-specific weighted sums;
- rank aggregation or Borda counts;
- count of top-\(N\) appearances across scenarios;
- Pareto-front selection;
- minimum-regret selection across scenarios;
- Bayesian model averaging when defensible priors are available.

For an early archival program, **top-\(N\) recurrence plus Pareto-front inclusion** is transparent and robust.

---

## 15. Recommended target-basket construction

A practical high-recall program should use partially overlapping baskets.

| Basket | Selection principle | Purpose |
|---|---|---|
| Distance core | Closest 20–30 components with minimal exclusions | Protect against over-modeling |
| Engineering favorable | Best distance-adjusted isolated, slowly rotating, dynamically quiet hosts | Represents a selective economical network |
| Natural graph neighbors | Recurring Delaunay, mutual-nearest, Gabriel, relative-neighborhood, or survivability-augmented Sun edges | Captures robust geometry |
| Selective hubs and gateways | Farther stars promoted after poor nearby hosts are filtered or because they open distinct regions | Covers backbone or picky-network architectures |
| Historical and kinematic | Past close approaches and persistent geometric neighbors | Covers old construction graphs |
| Exotic architecture | Isolated white dwarfs or other unusual high-upside lenses | Preserves model diversity |
| Pair-driven additions | Stars frequently appearing in strong redundancy pairs | Accounts for joint topology |

A reasonable initial portfolio might contain roughly 40–70 component-specific focal tracks after deduplication, depending on archive cost and overlap among binary-component loci.

Each retained candidate should carry one or more inclusion labels, for example:

- `DISTANCE_CORE`
- `ROBUST_GEOMETRIC_NEIGHBOR`
- `SELECTIVE_NETWORK_NEIGHBOR`
- `DENSE_HUB`
- `GATEWAY_HUB`
- `HISTORICAL_PARENT_CANDIDATE`
- `ENGINEERING_FAVORABLE`
- `COMPACT_LENS_WILDCARD`
- `PAIR_COMPLEMENT`
- `ARCHIVE_EASY`

These labels make the target list auditable and help interpret null results.

---

## 16. Suggested data fields

### 16.1 System-level fields

- system identifier and aliases;
- barycentric position and velocity;
- distance and covariance;
- multiplicity architecture;
- system age estimate;
- number and distances of neighboring eligible systems;
- Delaunay/Voronoi and other graph memberships;
- hub, gateway, and shortcut metrics;
- historical-neighbor persistence;
- candidate-system edge score under each scenario.

### 16.2 Component-level fields

- component identifier;
- mass and radius estimates;
- spectral and luminosity class;
- rotation period or \(v\sin i\);
- activity indicators;
- companion-induced acceleration estimates;
- astrometric acceleration and data-quality flags;
- \(z_{\min}\) and operating-distance scenarios;
- lens-host score and uncertainty;
- component-selection probability within the system.

### 16.3 Archive-level fields

- antipodal focal-line coordinates by epoch;
- search tube as a function of assumed node distance;
- archive names and observation intervals;
- wavelength and signature class;
- sensitivity and cadence;
- confusion and masking metrics;
- localization coverage fraction;
- expected search cost;
- completed-search and result status.

---

## 17. Uncertainty and robustness

### 17.1 Monte Carlo catalog realizations

Sample astrometric and stellar-parameter uncertainties, then recompute rankings. Report:

- mean or median rank;
- probability of top-\(N\) inclusion;
- confidence interval on graph metrics;
- probability of surviving each engineering filter;
- component-selection uncertainty in multiple systems.

### 17.2 Missing-data treatment

Use explicit imputation policies:

- neutral prior with a wide uncertainty interval;
- population-based prior conditional on stellar class;
- conservative penalty only for archival localization;
- no bonus for “no known companion” without detection-completeness information.

### 17.3 Sensitivity to assumptions

At minimum, vary:

- network age;
- persistent versus rewired links;
- maximum economic link distance;
- strength of the close-binary penalty;
- fraction of stars considered eligible;
- operating-distance policy;
- communication wavelength;
- importance of direct stellar power;
- relative weight of latency versus construction cost.

A candidate is especially valuable when it remains highly ranked under broad variations.

---

## 18. Common mistakes to avoid

1. **Treating binary components as independent graph destinations.** Use one system node, then expand to component-specific lens tracks.
2. **Dropping the nearest stars too early.** Engineering penalties may be negligible to an advanced civilization or irrelevant to historical parent links.
3. **Using proper motion as a direct stationkeeping cost.** Prefer angular acceleration, full relative motion, and ephemeris predictability.
4. **Assuming the target’s mass or radius changes the Sun-side focal minimum.** It changes the remote endpoint, not the solar endpoint.
5. **Assuming minimum focal distance is the actual operating distance.** A civilization may standardize operations much farther out.
6. **Rewarding unknown systems as though they were measured to be clean.** Unknown companion status is uncertainty, not evidence.
7. **Using known planet count without correcting for detection bias.** Separate resource hypotheses from dynamical penalties.
8. **Equating minimum degree two with strong network redundancy.** Pair and cycle value matter even in heuristic rankings.
9. **Promoting a generic graph hub without a direct Sun-edge argument.** Require adjacency, distance, gateway, or shortcut justification.
10. **Combining physical plausibility and archive coverage too early.** Preserve both scores so a poor archive does not masquerade as an implausible edge.

---

## 19. Minimal pre-simulation workflow

1. **Build a system catalog.** Consolidate close multiple components into systems while retaining component properties.
2. **Create the distance core.** Preserve the nearest 20–30 component focal tracks.
3. **Calculate host-quality terms.** Include companion acceleration, compactness, rotation, activity, and evolutionary state.
4. **Construct static geometric graphs.** Use several nearest-neighbor and proximity-graph definitions.
5. **Repeat after host filtering.** Identify farther systems that become natural neighbors in a picky network.
6. **Calculate hub and gateway proxies.** Include local density, incremental reach, and shortcut value.
7. **Propagate stellar geometry.** Evaluate historical neighbor status and edge persistence over plausible network ages.
8. **Estimate arrival-direction share.** Produce parent-link candidates without running a full expansion simulation.
9. **Rank candidate pairs.** Add stars favored by cheap cycles or independent-branch redundancy.
10. **Build signature-specific archive scores.** Calculate focal-line tracks, coverage, sensitivity, and confusion.
11. **Select the union of scenario baskets.** Deduplicate at the system and track levels while retaining inclusion labels.
12. **Publish both conditional and aggregate ranks.** Do not hide architecture dependence behind one number.

---

## 20. Recommended default priorities

Before catalog-specific analysis, a defensible ordering of criteria is:

1. **Current and plausible construction-era distance.**
2. **Natural geometric adjacency to the Sun.**
3. **Close-companion reflex acceleration and unresolved multiplicity risk.**
4. **Filtered adjacency after unfavorable nearby systems are removed.**
5. **Kinematic persistence under full three-dimensional relative motion.**
6. **Hub, gateway, shortcut, and pair-complement value.**
7. **Focal-distance and operating-environment engineering.**
8. **Rotation and activity, with wavelength-dependent weights.**
9. **Evolutionary stability, mainly as an extreme-object filter.**
10. **Resource indicators and known planetary architecture as weak scenario-specific terms.**
11. **Archive coverage and localization as a separate final multiplier.**

The expected outcome is not that the closest 20 stars are replaced. Rather, most remain in a protected core while a smaller number of more distant systems are promoted because they are unusually clean lens hosts, robust geometric neighbors, or plausible access points to the wider network.

---

## References

[^1]: S. G. Turyshev and V. T. Toth, work on the optics and focal region of the Solar Gravitational Lens; see, for example, *Resolved imaging of exoplanets with the Solar Gravitational Lens*, arXiv:2301.07495, and the associated SGL literature. <https://arxiv.org/abs/2301.07495>

[^2]: S. Kerby and J. T. Wright, *Stellar Gravitational Lens Engineering and SETI*, arXiv:2109.08657. The paper examines relay dynamics, companion-induced perturbations, rotation, and SETI implications. <https://arxiv.org/abs/2109.08657>

[^3]: N. Tusay et al., *A Search for Radio Technosignatures at the Solar Gravitational Lens Targeting Alpha Centauri*, *The Astronomical Journal* 164, 116 (2022), DOI: 10.3847/1538-3881/ac8358, arXiv:2206.14807. <https://arxiv.org/abs/2206.14807>

Additional background:

- M. Gillon, *A novel SETI strategy targeting the solar focal regions of the most nearby stars*, arXiv:1309.7586. <https://arxiv.org/abs/1309.7586>
- G. W. Marcy, S. K. Tellis, and E. H. Wishnow, *Laser Communication with Proxima and Alpha Centauri using the Solar Gravitational Lens*, arXiv:2110.10247. <https://arxiv.org/abs/2110.10247>
- Astropy documentation, *Accounting for Space Motion*, for propagating catalog coordinates under a constant-space-motion approximation. <https://docs.astropy.org/en/latest/coordinates/apply_space_motion.html>

---

## Status

This framework is intended for **pre-simulation target selection**. It is deliberately heuristic. Dynamic network simulations should later test which criteria genuinely predict Sun links under decentralized expansion, delayed information, redundancy construction, failures, and long-term rewiring.
