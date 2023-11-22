# Contributing to scpcaTools

## Coding style

We try to follow the [tidyverse style](http://style.tidyverse.org/).
You can use the [`styler` package](https://styler.r-lib.org) to perform styling locally, or use the pre-commit hooks described below.

If code needs styling during/after a PR, one of the repository owners can perform styling automatically by making a comment with just the word `/style` in the PR.
(Note that for this to work, the commenter must have set their membership to "Public" in https://github.com/orgs/AlexsLemonade/people)

## Pre-commit hooks

For convenience, we have included a set of [pre-commit hooks](https://pre-commit.com/) that can be used to automatically format code according to the above specifications, as well as to spell check and check for other common errors.

To use these hooks, install the `pre-commit` package according to your favorite method (`pip install pre-commit` or `conda install pre-commit`), then run `pre-commit install` in the `scpcaTools` directory.
This will install the hooks in the `.git/hooks` directory, and they will be run automatically when you commit changes.
If any of the hooks fail, the commit will be aborted, and you will need to fix the errors and re-commit.

Notably, the `spell-check` hook will report spelling errors, but will also add any words it finds to the dictionary file (stored at `inst/WORDLIST`).
This is convenient for many cases (where the word is real but unknown), but be sure to remove truly misspelled words from the dictionary file before committing, or they will not be caught in the future!
