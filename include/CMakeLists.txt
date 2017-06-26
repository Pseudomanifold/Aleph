# Install rule for all headers shipped by Aleph. This does not yet
# include the generated ones, though.
INSTALL(
  DIRECTORY   aleph
  DESTINATION include
  COMPONENT   library
  FILES_MATCHING PATTERN "*.hh"
)
