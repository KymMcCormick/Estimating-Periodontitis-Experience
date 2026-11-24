# Contributing

This repository is currently under active development.

## Style

- Use `snake_case` for derived variables and functions.
- Preserve original NHANES source variable names until the point where derived variables are created.
- Use explicit model iteration names rather than `new`, `final`, or `test`.
- Keep generated data and figures out of Git.

## Before changing model assumptions

Document any change to:

- extraction attribution schedules,
- cumulative attribution logic,
- missing-tooth burden assignment,
- tooth inclusion / exclusion rules,
- survey weighting, or
- diagnostic validation outcomes.

Add the change to `docs/model-specification-iteration01.md` or create a new iteration-specific specification file.
