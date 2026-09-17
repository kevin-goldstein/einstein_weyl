# Upload instructions

The bundle is ready to publish without the manuscript source or PDF.
Choose either GitHub or Zenodo as the public repository.
If you want a readable source repository and a preserved archival version, use GitHub first and then connect the GitHub release to Zenodo.

Before uploading, choose a license and replace the placeholder repository URL in `CITATION.cff`.
No license has been chosen on the authors' behalf.

## Option 1: publish on GitHub

1. Extract `einstein_weyl_public_repository_bundle.zip` into a new local directory.
2. On GitHub, create an empty **public** repository.
   Do not initialize it with a README, license, or `.gitignore`, because the bundle already contains the repository files.
3. In a terminal, enter the extracted directory and run:

   ```sh
   git init -b main
   git add .
   git commit -m "Publish computational certificate"
   git remote add origin https://github.com/OWNER/REPOSITORY.git
   git push -u origin main
   ```

   Replace `OWNER/REPOSITORY` with the repository you created.
   GitHub also supports `gh repo create --source=. --public --remote=origin --push` after the local commit.
4. Open the public repository in a signed-out browser window and confirm that `README.md`, `CITATION.cff`, `ancillary/`, and `fragments/` are visible.
5. Replace the manuscript placeholder

   ```tex
   \newcommand{\PublicRepositoryURL}{https://github.com/USERNAME/einstein-weyl-black-hole-proof}
   ```

   with the actual repository URL, rebuild the manuscript, and check the link in Appendix E.

These steps follow GitHub's official guides for [creating a repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository) and [adding locally hosted code](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github).

## Option 2: publish directly on Zenodo

1. Sign in to Zenodo, select **New upload**, and upload `einstein_weyl_public_repository_bundle.zip` as a single file.
   Zenodo recommends a ZIP archive for deposits with more than 20 files.
2. Use resource type **Software** and the title **Computational certificate for an Einstein--Weyl black hole**.
3. Add the creators in this order:
   - Kevin Goldstein
   - Vishnu Jejjala
4. Give both creators these affiliations:
   - Mandelstam Institute for Theoretical Physics, School of Physics, University of the Witwatersrand, Johannesburg 2050, South Africa
   - National Institute for Theoretical and Computational Sciences, South Africa
5. Use a short description such as: “Production sources, exact inputs, tests, verification scripts, successful logs, and generated numerical data supporting a computer assisted existence proof for a non-Schwarzschild black hole in Einstein--Weyl gravity.”
6. Add the keywords `Einstein--Weyl gravity`, `black holes`, `computer assisted proof`, `interval arithmetic`, and `Poincaré--Miranda theorem`.
7. Select the authors' chosen license, set the files to public, save the draft, and preview the record.
8. Publish the record after checking the file and metadata.
   Zenodo registers a DOI on publication, but the manuscript may simply link to the public record landing page.
9. Replace `\PublicRepositoryURL` in the manuscript with that landing-page URL, rebuild the manuscript, and check the link in Appendix E.

Zenodo's official instructions cover [creating an upload](https://help.zenodo.org/docs/deposit/create-new-upload/), [file preparation](https://help.zenodo.org/docs/deposit/manage-files/), and [creator metadata](https://help.zenodo.org/docs/deposit/describe-records/creators/).

## Optional: archive a GitHub release in Zenodo

After publishing the GitHub repository, connect the GitHub account in Zenodo, enable the repository, and create a GitHub release.
Zenodo then archives enabled repository releases.
The included `CITATION.cff` provides author and title metadata; update its repository URL before creating the release.

See Zenodo's official guides to [enabling a GitHub repository](https://help.zenodo.org/docs/github/enable-repository/) and [describing software](https://help.zenodo.org/docs/github/describe-software/).
