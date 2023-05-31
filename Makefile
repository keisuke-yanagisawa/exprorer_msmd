unittest:
	coverage run -m pytest . || echo ""

coverage: unittest
	coverage html