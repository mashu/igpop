# Releasing

The repository has three GitHub Actions workflows.

| Workflow | Trigger | What it does |
|---|---|---|
| `ci.yml` | push to `main`, pull requests | `cargo fmt --check`, `cargo clippy -D warnings`, `cargo test`. Tests do not use the network. |
| `release.yml` | pushing a tag `v*` | Builds release binaries for Linux x86_64, macOS arm64 and x86_64, and Windows x86_64, and attaches them to a GitHub release with generated notes. |
| `docs.yml` | push to `main` touching `docs/` | Builds this book with mdBook and publishes it to GitHub Pages. |

## One-time setup

In the repository settings, under *Pages*, set *Source* to
**GitHub Actions**. The site then appears at
<https://mashu.github.io/igpop/>.

## Cutting a release

```bash
# 1. set the version in Cargo.toml, commit
git tag v0.10.0
git push origin main v0.10.0
```

## Building the docs locally

```bash
cargo install mdbook
mdbook serve docs      # http://localhost:3000
```
