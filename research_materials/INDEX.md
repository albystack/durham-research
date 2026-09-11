# Research materials

This area is for the material surrounding the computational project. Its raw
contents are private/local by default so that publishing the code repository
does not accidentally publish personal correspondence or copyrighted papers.

```text
research_materials/
├── emails/             exported supervisor correspondence (`.eml`, `.pdf`, etc.)
├── supervisor_notes/   handwritten notes and diagrams supplied for the project
├── papers/             local reading copies; prefer links/citations in Git
└── bibliography/       BibTeX, citation notes, and a reading index
```

Recommended practice:

1. Keep raw email exports and private annotations local.
2. Record a short dated summary of each scientific decision in `docs/`.
3. Store bibliographic metadata or stable links rather than redistributing a
   paper unless its license permits that.
4. If a supervisor-provided note should become public, obtain permission and
   stage that file deliberately rather than overriding the ignore rules in
   bulk.

Professor Chhita's square-grid model instructions are already implemented as
an auditable mathematical contract in
[`../square_glauber_python/PERFECT_SAMPLING.md`](../square_glauber_python/PERFECT_SAMPLING.md)
and in the height/kernel source tests. The raw correspondence remains useful
context but is not the authoritative executable specification.
