For developers who want to contribute to the website/documentation of libigl.
If you want to preview changes to the libigl website before a commit, you can follow the instructions below.

1. Install mkdocs and the material theme
   ```bash
   pip install -U --user mkdocs mkdocs-material
   ```
2. Preview the website locally (in the root folder of the libigl project):
  ```bash
  python -m mkdocs serve
  ```
3. Build the website to generate the html locally (optional):
  ```bash
  python -m mkdocs build
  ```
4. Deploy the website directly to github (will overwrite the gh-pages branch of the remote repository):
   ```
   python -m mkdocs gh-deploy
   ```

!!! warning
    The `gh-deploy` script will overwrite the content of the `gh-pages` in the remote repository. Be sure of what you are doing before pushing new content with this command.

!!! tip
    Be careful to not have any `<>` characters in your email in your `.gitconfig`, otherwise the `gh-deploy` script will fail.

## References

- [MkDocs](http://www.mkdocs.org/)
- [Material Theme](https://squidfunk.github.io/mkdocs-material/)
