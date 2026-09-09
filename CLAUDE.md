# pointgroup-tables

## Agent skills

### Issue tracker

Issues live as GitHub issues in `hydeik/pointgroup-tables`, driven by the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

The five canonical triage roles, used verbatim as label strings. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.


## Language Instructions

The user prefers Japanese conversation.

Respond to the user in Japanese. For specialized technical terms, include the English term in parentheses where useful.

All generated or modified project artifacts must be written in English. This includes source code, comments, docstrings, tests, README files, documentation, configuration comments, commit messages, and pull request descriptions.

Do not include Japanese text in generated files unless explicitly requested by the user.

When explaining code changes, use Japanese in the chat, but keep code blocks and file contents in English.
