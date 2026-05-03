Notes to maintainers for building releases.
This is best done as a release PR so we can check it first.

1. Places to update the version number include:

   - htscodecs/htscodecs.h (used for program introspection)

   - configure.ac AC_INIT macro

   - configure.ac VERS_CURRENT, VERS_REVISION and VERS_AGE variables.
     See the long comment above for instructions of how these change.

   - NEWS files.


2. Ensure NEWS and README files are up to date.  NEWS is a git log
   summary.  README likely doesn't change unless something major needs
   mentioning.

   - At time of merging, set the date at the top of NEWS.


3. Test it all.
   - Push to github PR so the CI can validate for us.

   - make distcheck
     This also makes the tarball htscodecs-${vers}.tar.gz.


4. Merge into master


5. Add an annotated tag with minimal message, eg:

   - git tag -a v1.1 -m v1.1


6. Push master and --tags upstream to github


7. Make a new release on github.

   - Title: "htscodecs ${vers}"

   - Message: this is just a copy of NEWS.
     It's already in Markdown format, but double check the preview panel.

   - Upload the tarball produced from distcheck to the assets.


8. Finally, consider updating any packages that use this as a
   submodule to ensure they have the latest tagged release.

   This will invariably help OS distributions keep their package
   dependencies neatly in sync.
