---
name: Release checklist
about: Prepare for a new release version of scpcaTools
title: Prepare for scpcaTools release vX.X.X
labels: release
assignees: ''

---

## Steps for a new release of `scpcaTools`

- [ ] Are all dependencies/issues expected for this release are resolved? Add any that are not resolved as dependencies for this issue.
- [ ] Did all automated tests and the Docker build in this repo pass in CI? (You may have to wait for the docker image to build. Check [Github actions](https://github.com/AlexsLemonade/scpcaTools/actions/workflows/build-docker.yaml) before proceeding to the next step
- [ ] Test the new docker build in `scpca-nf`
  - [ ] Create a branch in [scpca-nf](https://github.com/AlexsLemonade/scpca-nf/) for testing the latest Docker version from the `main` branch.
  - [ ] Update [containers.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/config/containers.config) with the `edge` version of the scpca-tools container: ```SCPCATOOLS_CONTAINER = 'ghcr.io/alexslemonade/scpca-tools:edge'``` and push the branch to github.
  - [ ] If necessary, update the `scpca-nf` workflow to accommodate changes in the `scpcaTools` package  
  - [ ] Perform a test run of the workflow with:
```nextflow run AlexsLemonade/scpca-nf -r "<TESTBRANCH>" -profile batch```
You may want to specify libraries or runs for testing.
