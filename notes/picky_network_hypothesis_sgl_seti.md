# The Picky Network Hypothesis
## A selective interstellar SGL network and implications for SETI target ranking

## Summary

The **Picky Network Hypothesis** assumes that an advanced civilization building a self-expanding interstellar communication network does **not** attempt to link every stellar system. Instead, it admits only systems that offer a compelling reason for inclusion.

A system is linked because at least one of the following is true:

1. **It is an unusually good engineering or resource site** — for example, a dynamically clean star with favorable stellar-gravitational-lens (SGL) properties and abundant accessible material.
2. **It has unusually high graph value** — for example, it is a gateway, bridge, redundancy anchor, or efficient hub connecting otherwise poorly connected regions.
3. **It is intrinsically interesting** — for example, it hosts compelling biosignatures, a technological civilization, or an unusual planetary or astrophysical system.

Under this hypothesis, the Solar System need not connect to the nearest stars. A nearby but unattractive close binary may be skipped, while a more distant quiet K or M dwarf, a useful gateway, or a biologically interesting system may receive a direct link.

This has an important consequence for SETI searches near the Sun's focal lines: **a nearest-stars-only target list is incomplete by construction**. A high-recall archival search should combine proximity targets with engineering-favorable stars, graph candidates, biosignature/science targets, and a small number of architecture-dependent wildcards such as isolated white dwarfs.

---

# 1. Conceptual model

Let each stellar system have three broad desirability scores:

\[
E_i = \text{engineering and resource value},
\]

\[
G_i = \text{graph value},
\]

\[
I_i = \text{scientific or biological interest}.
\]

The key assumption is that selection behaves more like an **OR rule** than a conventional weighted average:

\[
\text{admit system } i \text{ if } E_i \text{ or } G_i \text{ or } I_i \text{ is exceptional}.
\]

A conceptual ranking score could therefore resemble

\[
S_i = \max(E_i, G_i, I_i)
+ \epsilon(E_iG_i + E_iI_i + G_iI_i),
\]

where the interaction term rewards systems that serve several roles simultaneously.

This is preferable to a simple score such as

\[
S_i = E_i + G_i + I_i,
\]

because a truly selective network may admit a mediocre relay if it hosts a compelling biosphere, or admit an uninteresting star if it is a uniquely valuable graph gateway.

The probability that the Sun is directly linked to candidate system \(i\) should then be treated separately:

\[
P(\odot \leftrightarrow i)
\propto
S_i\,A_{\odot i}\,R_{\odot i},
\]

where

- \(A_{\odot i}\) represents the geometric/economic convenience of attaching the Solar System to that node;
- \(R_{\odot i}\) represents the redundancy or routing value of the edge.

Thus, a highly desirable star still may not be a direct Solar neighbor if it can be reached more efficiently through another node.

---

# 2. System-level and component-level ranking

The Picky Network Hypothesis requires distinguishing three related but different questions.

## 2.1 System-level inclusion

Would the civilization choose to industrialize and connect the stellar **system** at all?

This depends on resources, graph position, planetary interest, long-term stability, and overall engineering value.

## 2.2 Component-level lens choice

If the system is multiple, which stellar component would be used as the gravitational lens?

Close binaries should not be treated as multiple independent network destinations merely because they contain multiple cataloged stars. The system is the graph node; individual stellar components are candidate lenses.

## 2.3 Archive-search priority

Even after estimating the probability of a physical link, archival detectability depends on additional factors:

- sky coverage;
- angular resolution;
- sensitivity;
- time baseline;
- astrometric uncertainty;
- confusion and background;
- expected station-distance range;
- signal type.

Therefore, **network-likelihood ranking and archival-searchability ranking should remain separate** until the final prioritization step.

---

# 3. Engineering and SGL-host value

## 3.1 Minimum focal distance

For a lens star of mass \(M_\star\) and radius \(R_\star\), the beginning of the gravitational focal line scales approximately as

\[
z_{\min} \propto \frac{R_\star^2}{M_\star}.
\]

A convenient compactness proxy is therefore

\[
Q_{\rm focal} \propto \frac{M_\star}{R_\star^2}.
\]

Lower \(R^2/M\) generally means the usable focal region begins closer to the star.

Important qualification: the focal line continues outward indefinitely. A civilization is not forced to operate exactly at \(z_{\min}\). Therefore, a very small minimum focal distance is primarily an **option value**, whereas a very large minimum distance is a genuine logistical burden.

For a Sun-other-star link, the candidate star's properties affect the endpoint near the candidate star. The Sun-side endpoint always uses the Sun as the lens, so its minimum focal distance is determined by the Sun, not by the remote target.

## 3.2 Multiplicity and reflex acceleration

A binary flag is too crude. The physically important quantity is the acceleration of the chosen lens star caused by companions:

\[
a_{\rm reflex}
\sim
\sum_k \frac{GM_k}{a_k^2},
\]

with orbital phase, eccentricity, inclination, and projection included where known.

Strongly unfavorable cases include:

- close stellar binaries;
- contact binaries;
- compact hierarchical multiples;
- stars with massive close companions;
- systems with large measured astrometric acceleration.

More favorable cases include:

- apparently single stars;
- extremely wide binaries;
- systems whose companions induce only very small accelerations at the selected lens component.

Close companions are one of the strongest astrophysically motivated reasons for a picky network to skip a nearby system.

## 3.3 Rotation and oblateness

Rapid stellar rotation can make the gravitational lens less axisymmetric. A selective civilization might prefer

- slow rotation;
- low projected rotational velocity;
- low oblateness;
- stable stellar structure.

This is probably a secondary factor for an extremely advanced civilization, but it is a plausible tie-breaker among otherwise similar candidates.

## 3.4 Activity and plasma environment

Magnetic activity, stellar wind, and coronal plasma may complicate high-gain communication, especially at radio wavelengths.

Potential ranking variables include

- X-ray luminosity;
- UV activity;
- H-alpha emission;
- flare rate;
- rotation period;
- age;
- wind proxies.

This penalty should be **wavelength-dependent**. An active M dwarf may be undesirable for a radio relay but less problematic for optical or near-infrared operation.

## 3.5 Evolutionary stability

For permanent or very long-lived network infrastructure, favorable systems include stars with long remaining stable lifetimes and slow structural evolution.

Broadly favorable classes:

- old, quiet K dwarfs;
- old, inactive M dwarfs;
- ordinary stable G dwarfs;
- isolated compact remnants under some architectures.

Broadly unfavorable classes:

- giants and subgiants;
- massive short-lived stars;
- accreting compact objects;
- interacting binaries;
- young, rapidly evolving systems.

Among ordinary FGKM dwarfs, lifetime may be less important than multiplicity, activity, and resources unless the expected network lifetime is extremely long.

---

# 4. Resource value

A relay is not just a lens; it is an industrial site. A picky civilization may favor systems where infrastructure can be constructed, fueled, repaired, and expanded cheaply.

Potentially desirable properties include

- large inventories of rocky and metallic material;
- asteroid or planetesimal belts;
- icy bodies rich in volatiles;
- giant planets containing large hydrogen inventories;
- substantial total planetary-system mass;
- debris belts or remnant planetesimal populations;
- sufficient stellar luminosity for inner-system power generation.

There is a tradeoff between resources and dynamical cleanliness. A giant planet may provide enormous useful mass while simultaneously inducing reflex acceleration on the lens star.

Thus, resource value and lens-dynamics value should be kept as separate terms.

Observationally, resource ranking is uncertain because nearby-system inventories are incomplete. Known planet count, known debris disks, and metallicity should therefore be treated as **weak or moderate proxies**, not direct measurements of usable material.

---

# 5. Expected stellar classes under the hypothesis

## 5.1 Old, quiet K dwarfs

These may be the best general-purpose candidates.

Advantages:

- long remaining lifetime;
- relatively compact compared with solar-type stars;
- moderate luminosity and abundant stellar power;
- typically quieter than young M dwarfs;
- potentially favorable planetary systems;
- good compromise between lens compactness and industrial convenience.

A picky backbone architecture may favor these disproportionately.

## 5.2 Old, inactive early- and mid-M dwarfs

These may be excellent if local starlight at the relay is not a limiting factor.

Advantages:

- small radii and favorable focal-distance scaling;
- extremely long lifetimes;
- common in the local neighborhood;
- many host terrestrial planets.

Disadvantages:

- low luminosity;
- activity can remain significant for long periods;
- flare and plasma environment varies strongly from star to star.

Thus, **old and inactive** is more important than merely being an M dwarf.

## 5.3 Ordinary G dwarfs

These are reasonable multipurpose nodes, especially when they host rich or biologically interesting planetary systems.

They are less optimized for minimum focal distance than smaller stars, but may offer

- abundant stellar power;
- mature planetary systems;
- familiar long-term stability;
- rich resource environments.

## 5.4 Isolated white dwarfs

These are important architecture-dependent wildcards.

Advantages:

- extremely small radius for substantial mass;
- very favorable formal focal-distance scaling;
- low optical luminosity;
- long-lived remnant state;
- potentially clean lens geometry if isolated and slowly rotating.

Disadvantages:

- little natural stellar power at large focal distances;
- post-main-sequence system history;
- uncertain resource availability;
- stronger gravity close to the remnant;
- accreting or binary white dwarfs can be highly unfavorable.

A high-recall SGL SETI program should explicitly include an **isolated-white-dwarf basket** rather than assuming ordinary main-sequence architectures only.

---

# 6. Graph value in a selective network

The graph should be constructed among **desirable systems**, not all stars.

A star surrounded by many nearby but unattractive binaries may have low value, whereas a somewhat farther star connecting several high-quality systems may become a major hub.

## 6.1 Dense hubs

A dense hub lies near many other desirable systems.

A simple metric is

\[
H_{\rm dense}(i)
=
\sum_j q_j f(d_{ij}),
\]

where \(q_j\) represents the desirability of system \(j\).

A good hub may support many short outgoing edges from one industrial base.

## 6.2 Gateway nodes

A gateway occupies a strategically important position between otherwise weakly connected regions.

Possible proxies include

- betweenness centrality in a filtered geometric graph;
- bridge-removal value;
- Delaunay adjacency among desirable systems;
- access to systems on the far side of a local void;
- reduction in total path length when included;
- reduction in required long edges.

This is the main mechanism by which a star 15–30+ light-years away could become a more plausible Solar neighbor than several stars inside 10 light-years.

## 6.3 Redundancy value

If each admitted system must have at least two links, some stars become valuable because they close important cycles or provide independent paths.

A direct Sun-star edge may be favored when it

- eliminates a bridge;
- creates a second route to a backbone region;
- connects the Solar System to a different spatial sector;
- reduces failure vulnerability;
- materially shortens alternate routes.

## 6.4 Kinematic persistence

An ancient network may favor systems whose useful graph relationships persist for long periods.

Prefer systems with

- low relative velocity with respect to neighboring nodes;
- slowly changing edge lengths;
- long-lived geometric adjacency;
- stable membership in the same local stellar neighborhood.

Raw proper motion alone is not sufficient. Full 3D relative motion and time-dependent neighborhood geometry are more meaningful.

---

# 7. Scientific and biological interest

A sufficiently interesting system may be admitted despite mediocre engineering properties.

## 7.1 Highest-interest systems

Potentially overriding all engineering penalties:

- confirmed technological civilization;
- strong technosignatures;
- unambiguous complex biosphere;
- robust atmospheric disequilibrium indicating biology;
- industrial atmospheric chemistry;
- artificial structures or anomalous energy use.

## 7.2 Strong science targets

Likely to receive high priority:

- temperate rocky planets with atmospheres;
- surface liquid-water evidence;
- several independent biosignature candidates;
- multiple potentially habitable planets in one system;
- unusual planetary architectures.

## 7.3 Moderate-interest science nodes

Possible reasons for inclusion even without life:

- ocean worlds;
- young terrestrial planets;
- resonant planetary systems;
- rare compact objects;
- unusual atmospheric chemistry;
- unusual debris architectures;
- systems useful for long-baseline astrophysical monitoring.

---

# 8. The redundancy penalty for science-only nodes

A key implication of a two-link minimum is that an interesting system is not cheap to add if it sits far from the existing network.

If an interesting system \(X\) must connect to two existing nodes \(A\) and \(B\), then approximately

\[
C_{\rm add}(X)
\approx
C_{AX}+C_{BX}.
\]

Thus, two equally interesting biosignature systems may have very different inclusion probabilities if one lies close to the backbone and the other requires two long edges.

This means scientific-interest admission is more realistically approximated by

\[
P(\text{science node})
\propto
I_i \times \text{attachment convenience}.
\]

An interesting system adjacent to two high-quality backbone stars may be much more likely to be connected than an equally interesting isolated system.

This also has a direct implication for the Solar System: if Earth is the reason the Solar System was admitted, the Sun's neighbors may be the nearest **pre-existing backbone nodes**, not the nearest stars overall.

---

# 9. Nearby examples that fit the hypothesis

A qualitative screen of roughly the nearest 100 stellar systems reveals several particularly interesting categories.

## 9.1 Strong multipurpose M-dwarf candidates

### GJ 887 / Lacaille 9352

Why it stands out:

- very nearby;
- relatively quiet M dwarf;
- compact SGL-host geometry;
- rich planetary system;
- habitable-zone interest.

This is a strong example of a system where engineering and science motivations reinforce one another.

### Ross 128

Why it stands out:

- old, relatively quiet M dwarf;
- compact lens host;
- terrestrial planet of high scientific interest;
- dynamically simpler than many nearby active multiples.

### GJ 1061

Why it stands out:

- low-mass star with favorable focal-distance scaling;
- multiple low-mass planets;
- relatively low activity compared with more troublesome flare stars.

### Teegarden's Star

Why it stands out:

- very low stellar mass and radius;
- relatively low activity for its class;
- multiple terrestrial planets;
- strong combined engineering/science case.

### Luyten's Star / GJ 273

Why it stands out:

- quiet nearby M dwarf;
- confirmed low-mass planet in or near the habitable zone;
- plausible combination of compact lens geometry, scientific interest, and minor-body resources.

## 9.2 Strong K/G candidates

### Epsilon Indi A

Why it stands out:

- K dwarf;
- long-lived stable host;
- known giant planet;
- substantial planetary/substellar architecture;
- likely strong resource value.

The giant planet introduces reflex acceleration, so this may be more attractive as an industrial/science node than as a perfectly clean lens.

### Tau Ceti

Why it stands out:

- mature nearby G star;
- extensive debris architecture;
- potential resource value;
- plausible planetary interest.

### Sigma Draconis

Why it stands out:

- mature nearby K dwarf;
- good candidate for a pure backbone role;
- not dependent on known spectacular planets to justify inclusion.

Its importance would rise substantially if filtered graph analysis shows high gateway or redundancy value.

### HD 20794 / 82 Eridani

Why it stands out:

- nearby G dwarf;
- multiple confirmed planets;
- one low-mass planet intersects or approaches the habitable-zone regime;
- strong science-target motivation.

This is exactly the kind of system that a distance-only nearest-star search can miss.

## 9.3 Compact-remnant wildcards

### van Maanen's Star

Why it stands out:

- nearby isolated white dwarf;
- extreme compact-lens geometry;
- evidence of planetary-material accretion;
- no close luminous stellar companion.

### LAWD 37 / GJ 440

Why it stands out:

- nearby isolated white dwarf;
- clean compact-remnant test case;
- useful for probing whether an SGL network strongly favors compact lenses.

These two systems deserve explicit inclusion in high-recall archival searches even if the assumed probability of compact-remnant relays is modest.

## 9.4 Science-overrides-engineering example

### GJ 876

Why it stands out:

- dynamically rich multiplanet system;
- famous resonant architecture;
- high scientific value.

Massive planets may make it a less pristine lens host, but the planetary system itself could justify inclusion.

---

# 10. Systems a picky network may skip

If a nearby system has **no exceptional science target and no major graph value**, a picky network may strongly downrank it for

- close binarity;
- compact hierarchical multiplicity;
- rapid rotation;
- strong astrometric acceleration;
- giant or subgiant evolutionary state;
- youth and strong magnetic activity;
- high relative velocity and transient neighborhood geometry;
- poor local connection to other desirable nodes;
- lack of accessible resources.

This creates a qualitatively different topology from a nearest-neighbor network.

A selective architecture can plausibly skip several stars inside 5–10 light-years while connecting directly to better systems in the 10–20+ light-year shell.

---

# 11. Recommendations for SGL SETI searches

The goal of a pre-simulation search program should be **high recall across competing selective-network architectures**.

Do not rely on a single ranked list.

## 11.1 Search basket A: nearest inclusive targets

Include the nearest approximately 20–30 stellar components regardless of engineering quality.

Purpose:

- protect against overfitting speculative assumptions;
- cover persistent parent links from an earlier expansion phase;
- capture civilizations that regard binary/rotation penalties as trivial.

## 11.2 Search basket B: engineering-favorable backbone candidates

Promote stars with combinations of

- low \(R^2/M\);
- single-star or dynamically quiet architecture;
- low companion-induced acceleration;
- slow rotation;
- low activity;
- long stable lifetime;
- substantial resource inventory.

Priority stellar classes:

- old quiet K dwarfs;
- old inactive early/mid-M dwarfs;
- selected stable G dwarfs.

## 11.3 Search basket C: graph/gateway candidates

Before full expansion simulation, construct static filtered graphs among desirable systems and promote stars that are

- Delaunay neighbors of the Sun;
- mutual nearest neighbors after poor hosts are removed;
- strong gateway nodes;
- bridge-removal or cycle-closing candidates;
- unusually persistent neighbors under stellar motion;
- high-value shortcuts between desirable regions.

This basket is the main defense against missing farther backbone nodes.

## 11.4 Search basket D: science-interest systems

Include nearby systems with

- temperate terrestrial planets;
- multiple low-mass planets;
- strong biosignature potential;
- unusual resonant architectures;
- unusual debris or planetary systems;
- compelling direct-imaging prospects.

These targets test the possibility that network membership is driven by scientific interest rather than relay engineering.

## 11.5 Search basket E: compact-remnant wildcards

Explicitly include nearby isolated white dwarfs such as

- van Maanen's Star;
- LAWD 37 / GJ 440.

These probe a very different SGL engineering architecture that conventional habitable-star target lists would miss.

---

# 12. Provisional high-priority target set

Before dynamic network simulation, the following systems deserve elevated attention under the Picky Network Hypothesis.

## Very high priority

- GJ 887 / Lacaille 9352
- Ross 128
- GJ 1061
- Teegarden's Star
- Luyten's Star / GJ 273
- Epsilon Indi A

These combine multiple independent reasons for inclusion.

## High priority

- Barnard's Star
- Lalande 21185
- Epsilon Eridani
- Tau Ceti
- HD 20794 / 82 Eridani

## Backbone candidates requiring graph analysis

- Sigma Draconis
- additional quiet single K dwarfs in the 10–25 light-year volume

These may rise sharply if filtered graph geometry gives them strong gateway or redundancy value.

## Architecture wildcards

- van Maanen's Star
- LAWD 37 / GJ 440

## Science-overrides-engineering

- GJ 876
- other unusually rich or dynamically interesting planetary systems

---

# 13. Recommended archival-search strategy

Because archival focal-line searches can be relatively cheap compared with building a full network simulation, the preferred strategy is to search the **union** of several independently motivated baskets.

A practical program could contain roughly 40–70 component-specific solar focal-line tracks:

- 20–30 nearest components;
- 10–15 additional engineering-favorable stars;
- 10–15 graph or gateway candidates;
- 5–10 science-interest or compact-remnant wildcards;
- extra components where a multiple system leaves ambiguity over which star would be used as the lens.

The exact number is less important than preserving diversity of hypotheses.

For each candidate, maintain separate quantities for

1. **network-edge prior** — probability that a Sun-star link exists;
2. **component-lens prior** — probability that this stellar component is used;
3. **focal-line localization** — predicted sky track as a function of station distance and epoch;
4. **archive coverage** — sensitivity, cadence, confusion, and temporal completeness;
5. **signature model** — thermal emission, optical/IR point source, radio emission, laser leakage, propulsion, or other technosignature.

Do not use archival coverage to modify the physical network prior. Combine the two only when deciding search order.

---

# 14. Recommended pre-simulation analyses

Before building a full autonomous expansion simulation, several inexpensive analyses can substantially refine the target list.

## 14.1 Engineering screen of the nearest ~100 systems

Estimate or compile

- mass;
- radius;
- \(R^2/M\);
- multiplicity;
- companion-induced acceleration;
- rotation;
- activity;
- evolutionary state;
- known planets;
- debris architecture;
- resource proxies.

## 14.2 Filtered static graph analysis

Construct graphs only among favorable or interesting systems.

Useful graph families:

- mutual k-nearest-neighbor;
- distance-threshold;
- 3D Delaunay;
- Gabriel graph;
- relative-neighborhood graph;
- minimum-spanning tree plus redundancy augmentation.

Repeat after removing increasingly large fractions of poor engineering hosts.

Stars that repeatedly become direct Sun neighbors are strong pre-simulation candidates.

## 14.3 Historical geometry

Propagate nearby stellar positions through plausible construction epochs and measure

- past Sun-star distance;
- closest approach;
- duration as a natural graph neighbor;
- persistence of Delaunay or mutual-neighbor relations.

This is especially useful for ancient persistent networks.

## 14.4 Neighbor-pair analysis

Because the Sun should have at least two connections, rank **pairs** as well as individual stars.

Two useful pair classes:

- **cheap-cycle pairs**, where the Sun and two candidate nodes form a low-cost triangle;
- **independent-branch pairs**, where the two neighbors connect the Sun to distinct backbone regions.

A moderately ranked star may become a high-priority archival target if it occurs frequently in high-value pairs.

---

# 15. Main predictions of the Picky Network Hypothesis

If this hypothesis is correct, an SGL network around the Sun should show several qualitative features.

## Prediction 1: nearby stars are not sufficient

The most probable Solar neighbors need not be the absolute nearest stars.

## Prediction 2: close binaries are underrepresented

Unless scientifically valuable or graph-critical, dynamically troublesome multiple systems should appear less often than their raw spatial density suggests.

## Prediction 3: quiet K and M dwarfs are overrepresented

Old, quiet, dynamically clean low-mass stars should appear more often than under a purely distance-driven network.

## Prediction 4: some direct Sun links may extend to 10–30+ light-years

These would likely correspond to gateway, backbone, or high-interest systems.

## Prediction 5: isolated white dwarfs may appear as unusual hubs

If compact-lens efficiency is strongly valued and local stellar power is not essential, isolated white dwarfs may be disproportionately useful.

## Prediction 6: known planetary-interest systems may be overrepresented

A selective exploratory network should preferentially include systems with terrestrial planets, biosignature potential, unusual planetary architectures, or other high scientific value.

## Prediction 7: the Solar System may itself be an interest node

Earth could be the principal reason the Solar System was admitted to the network. In that case, our direct neighbors may be nearby backbone nodes chosen to attach an interesting Solar System cheaply and redundantly.

This possibility is particularly important for SETI because it predicts that our SGL neighbors may be **good relays rather than Sun-like stars**.

---

# 16. Bottom line

The Picky Network Hypothesis replaces the assumption

> every nearby star is a potential network node

with

> only unusually useful, strategically positioned, or scientifically interesting systems merit infrastructure.

Under this model, the most promising Sun-neighbor candidates are not simply the nearest stars. The search should emphasize four overlapping populations:

1. **nearby unavoidable candidates**, retained for robustness;
2. **quiet, dynamically clean K/M backbone stars**;
3. **farther gateway nodes identified by filtered graph geometry**;
4. **scientifically exceptional systems and compact-remnant wildcards**.

For SGL SETI, the strongest practical recommendation is therefore to conduct broad archival searches across a diversified target portfolio before waiting for a full network-expansion simulation. The local 10–20+ light-year shell already contains several systems that appear more purposeful under a selective-network architecture than many of the absolute nearest stars.

