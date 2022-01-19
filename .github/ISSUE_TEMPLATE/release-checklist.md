---
name: Release checklist
about: Prepare for a new release version of scpcaTools
title: Prepare for scpcaTools release vX.X.X
labels: release
assignees: ''

---

## Steps for a new release of `scpcaTools`

### Preparing for the release and testing

- [ ] Are all of the issues planned for this release resolved? If there are any issues that are unresolved, mark this issue as blocked by those on ZenHub.
- [ ] Did all automated tests and the Docker build in this repo pass in CI? (You may have to wait for the docker image to build. Check [Github actions](https://github.com/AlexsLemonade/scpcaTools/actions/workflows/build-docker.yaml) before proceeding to the next step
- [ ] Test the new docker build in `scpca-nf`
  - [ ] Create a branch in [scpca-nf](https://github.com/AlexsLemonade/scpca-nf/) for testing the latest Docker version from the `main` branch.
  - [ ] Update [containers.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/config/containers.config) with the `edge` version of the scpca-tools container: ```SCPCATOOLS_CONTAINER = 'ghcr.io/alexslemonade/scpca-tools:edge'``` and push the branch to github.
  - [ ] If necessary, update the `scpca-nf` workflow to accommodate changes in the `scpcaTools` package  
  - [ ] Perform a test run of the workflow with:
```nextflow run AlexsLemonade/scpca-nf -r "<TESTBRANCH>" -profile batch```
If the [default `run_ids`](https://github.com/AlexsLemonade/scpca-nf/blob/main/main.nf#L4-L9) for the workflow do not cover expected changes in this tools package, you may want to specify particular test samples or projects with the `--run_id` option (use your judgement).

### Creating a release
- [ ] On the [releases page](https://github.com/AlexsLemonade/scpcatools/releases), choose `Draft a new release`.
- [ ] In `Choose a tag`, type a new release number using semantic versioning (vX.X.X) (you did update the title of this issue to match, right?), then click `Create a new tag: vX.X.X on publish`.
- [ ] Write a description of the major changes in this release. You may want to start with the Auto-generated release notes to save time.
- [ ] Save a draft to return to later if testing isn't done yet, otherwise:
- [ ] Publish the release!
