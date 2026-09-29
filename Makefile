.PHONY: test demo-report install-app app

PYTHON ?= python3

test:
	$(PYTHON) -m unittest discover -s tests -v

demo-report:
	$(PYTHON) scripts/build_competition_demo_report.py

install-app:
	./scripts/install_app.zsh

app:
	./scripts/run_app.zsh
