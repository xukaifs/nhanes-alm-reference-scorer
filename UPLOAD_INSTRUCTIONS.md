# Upload instructions for v2.2.0

This is an **update/merge package** for the existing GitHub repository, not a replacement repository snapshot.

1. Extract this ZIP locally.
2. Copy/merge its contents into your existing local clone of `xukaifs/nhanes-alm-reference-scorer`.
3. If prompted about the two existing `models/*.rds` files, **keep the existing GitHub versions**. This package intentionally does not include them.
4. Review the diff before committing. In particular, the two frozen RDS model files must remain present after the merge.
5. Commit the Paper 1B resources, push to `main`, then create the `v2.2.0` release.

Recommended commit message:

`Add Paper 1B cross-population validation resources for v2.2.0`
