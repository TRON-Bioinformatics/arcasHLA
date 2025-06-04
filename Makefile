all: format test

.PHONY: format lint test

# Requires dev dependencies to be installed.
format:
	# Markdown files.
	fd ".*\.md$$" -H --exec mdformat {}

	# Python scripts.
	black scripts/*.py

# Requires dev dependencies to be installed.
lint:
	# Markdown files.
	fd ".*\.md$$" -H --exec mdformat --check {}

	# Python scripts.
	black --check scripts/*.py

test:
	pytest
