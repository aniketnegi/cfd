# Engineering and documentation standard

## Aim

Write software that a competent engineer can understand, change, test, and review with limited context.

Optimize software for these qualities:

- Correctness
- Conceptual integrity
- Local reasoning
- Explicit state and dependencies
- Change locality
- Readable control flow
- Stable, narrow contracts

Do not optimize for line counts, file counts, function counts, class counts, or design-pattern counts.

These rules are defaults that support this aim. Do not apply them mechanically.
If an exception makes code clearer or safer, make the exception. State the trade-off briefly.

## Prepare a code change

1. Before you propose a design, inspect the relevant implementation, tests, public contracts, and repository conventions.
2. Identify these facts:
   - Behavior that will change
   - Invariants that must remain true
   - Owner of the affected data and policy
   - Applicable side effects
   - Compatibility constraints
   - Verification method
3. For nontrivial work, present a plan in five lines or fewer.
4. Make the smallest coherent change. Do not redesign, rename, reformat, or clean unrelated code.
5. Ask for clarification only if ambiguity can change behavior, destroy data, weaken security, or change a public contract.
6. Otherwise, inspect the repository and proceed.

## Modules and boundaries

A module has one coherent responsibility. A module is not necessarily a file, class, package, service, or directory.

Before you add a module, layer, interface, wrapper, factory, helper, or file, answer these questions:

1. What decision, invariant, policy, lifecycle, authority, external effect, or knowledge does it own?
2. What complexity does it hide?
3. What likely change does it keep local?
4. Why is its interface, call, or navigation cost useful?

If these questions have no clear answer, do not add the boundary.

Separate responsibilities by knowledge ownership and reasons to change. Do not separate them only by execution order.

An interface should be substantially simpler than the implementation that it hides.

Each layer must add policy, invariant enforcement, ownership, effect isolation, translation, lifecycle control, or a useful capability.
Collapse a layer that only forwards arguments and results.

Keep an interface when it creates a stable contract or isolates unstable, external, privileged, or difficult-to-test code.
Remove an interface that supports only hypothetical substitution.

Create abstractions from concrete evidence.
Put two code units behind one abstraction only when they meet all these conditions:

- They have the same meaning
- They obey the same invariant
- They change for the same reason

Occurrence count does not prove a shared abstraction. Temporary duplication is safer than a guessed abstraction.

Point dependencies toward stable domain contracts.
Keep databases, networks, filesystems, operating-system calls, user interfaces, and vendor APIs at explicit edges.
Avoid dependency cycles.

## Functions, files, and control flow

Do not use a fixed line limit for functions or files.

Keep a cohesive, sequential algorithm together when extraction causes one or more of these results:

- Obscures execution order
- Scatters shared state
- Enlarges the parameter surface
- Adds single-use calls
- Forces readers to move between files

Extract a part when it meets one or more of these conditions:

- Has a meaningful contract
- Isolates a side effect
- Reduces state or nesting depth
- Represents a domain operation
- Supports independent changes and tests

Use a small helper when its name establishes a useful domain concept or invariant.
Do not use a helper only to shorten its caller.

Split a file when it contains responsibilities that change for independent reasons.
Merge files when readers must navigate repeatedly and boundaries add no policy, invariant, or reusable capability.

Do not require one class or type in each file.
Keep small, related types and functions together when this structure makes code easier to read.

Prefer direct call paths. Identify the semantic work at each additional call.
Collapse calls that only rename or forward data.

Reduce deep nesting with guard clauses, better data structures, or explicit state machines when they improve clarity.
Do not change direct code only to satisfy a complexity metric.

## Data, state, and types

Design the data model, ownership, and invariants before you design control flow.

Make dependencies and mutable state explicit.
Avoid hidden globals, implicit registries, service locators, and ambient context unless the platform requires them.

Prefer pure transformations for computation. Keep mutation local.
Do not force functional style when controlled in-place mutation is clearer or necessary for measured performance.

Use types when they express domain meaning, encode ownership, or exclude invalid states.
Do not create decorative wrapper types that add no check, behavior, or clarity.

Validate untrusted data at system boundaries.
After a boundary establishes an invariant, internal code may rely on it without repeated defensive checks.

Do not swallow errors.
An error should identify the attempted operation, relevant object or value, and underlying cause.
Preserve useful context when you propagate an error.
Make sure that each failure path releases or rolls back resources that it owns.

Do not add performance complexity without measurements.
You can make an exception for a documented hard constraint.
When you optimize, report the benchmark and trade-off.

## Tests

Test observable behavior, contracts, and invariants. Do not test private implementation structure.

Add a regression test for each fixed defect when you can reproduce the defect reliably.

Test boundary conditions and important failure paths. Do not test only the expected path.

Keep tests deterministic unless randomness is the test subject.
Make external dependencies explicit and controlled.

Run focused tests first. Then run the relevant wider suite.
Run applicable formatters, linters, type checkers, and static analysis tools.
State exactly what you ran. Never claim success for a command that you did not run.

## Comments

Comments preserve information that code cannot express clearly.

Use comments for these subjects:

- Rationale and rejected alternatives
- Invariants, ownership, and lifecycle
- Preconditions, postconditions, and side effects
- Nonobvious compatibility or performance constraints
- Unusual error handling
- Mathematical, protocol, regulatory, or algorithmic sources
- Reasons that an apparently simpler implementation is incorrect

Do not comment visible mechanics, repeat a signature, or narrate each statement.
Do not use comments to compensate for unclear code. Rewrite the code first.

Document public contracts and surprising internal contracts.
Do not require the same comment template for every function.

Delete stale comments and commented-out code. Version control is the archive.

## Documentation

Before you create pages, plan the documentation set.
Identify audiences, reader tasks, required concepts, and information that readers must find.

Give each page one dominant reader need:

- Tutorial: Guide a reader to a first successful result.
- How-to guide: Help a reader complete a real task.
- Reference: Give exact behavior and contracts.
- Explanation: Explain design, rationale, limits, and trade-offs.

These categories organize documentation. They are not strict content rules.
Include enough local context to prevent unnecessary page changes.

A page must complete the task or answer the question in its title.
Use links for optional details, adjacent tasks, or full reference material.
Do not make readers follow a link to get an essential step.

Prefer a small number of substantial, navigable pages over many fragments.
Split a page when it serves independent reader needs or maintenance boundaries.
Do not split a page only because it has multiple headings.

Use cross-links deliberately:

- Tutorials lead readers to the next practical tasks.
- How-to guides lead readers to exact reference entries.
- Reference entries lead readers to related usage and rationale.
- Explanations lead readers to the procedures and contracts that they discuss.

Use one stable term for each concept.
State limitations, unsupported uses, compatibility constraints, and trade-offs directly.

Make examples correct, representative, and safe to copy.
Show failure behavior when it is important.

## Prose

Write direct technical prose.

Lead with the fact, decision, or action. Use concrete subjects and active verbs.
Keep subjects near their verbs. Give known context before new information.

Remove these forms:

- Ceremonial introductions
- Generic praise
- Marketing language
- Vague quality claims
- Repeated summaries
- False contrasts such as “not only X but also Y”
- Decorative lists of three
- Unnecessary headings
- Rhetorical questions that delay the answer
- Commentary about what the document will discuss
- Conclusions that repeat the prior section

Use approved STE words when they preserve technical meaning.
Treat necessary software and research terms as technical terms.
Do not replace a precise technical term with an inaccurate approved word.

Do not use adjectives such as `robust`, `scalable`, `simple`, or `flexible` without evidence.
State the applicable failure mode, load, operation, or design axis.

Do not imitate conversational LLM patterns.
Do not praise the question, announce an exploration, or narrate the writing process.

## Review a change

Review nontrivial work in separate passes.

### Review scope and architecture

- Confirm that the change is at the correct boundary.
- Confirm that each unit owns one coherent responsibility.
- Check for new dependency cycles and misplaced dependencies.

### Delete and flatten

- Remove unnecessary code, layers, options, wrappers, parameters, and configuration.
- Collapse calls that only forward data.
- Remove speculative scaffolding.

### Verify correctness

- Verify all applicable invariants.
- Check boundary cases, failure paths, security checks, ownership, and resource cleanup.
- Verify required compatibility.

### Check readability

- Use truthful and consistent names.
- Make state explicit.
- Make control flow visible.
- Remove unnecessary file navigation.

### Check tests and documentation

- Confirm that tests establish the intended behavior.
- Confirm that comments preserve rationale instead of mechanics.
- Update affected public contracts and user documents.
- Remove filler from prose.

### Read the final diff

- Read the final diff without implementation-process context.
- Remove content that reviewers do not need for the logical change.
- Verify each reported command and result.

## Numerical and research software

Treat numerical assumptions as part of the public contract.

Record these facts when they are relevant:

- Physical units and dimensions
- Coordinate systems and sign conventions
- Array shapes, data types, layouts, and valid ranges
- Precision assumptions
- Tolerances and convergence criteria
- Random seeds and nondeterminism sources
- Dataset, model, and parameter provenance
- Hardware-dependent or compiler-dependent behavior

Preserve operation order when reassociation can change stability or reproducibility.

Validate numerical work with the strongest available evidence:

- Analytical or manufactured solutions
- Trusted reference implementations
- Experimental measurements
- Conservation laws and physical invariants
- Dimensional checks
- Convergence and sensitivity studies
- Limiting and degenerate cases

Separate the mathematical model, discretization, solver, and input or output only when each owns a decision that can change.
Do not create modules from pipeline steps only because they run in sequence.

Keep measured computational kernels flat and aware of data layout when abstractions add material cost.
Record the applicable benchmark, precision, hardware, and input scale.

Put equation, paper, standard, or derivation references near the code that depends on them.
