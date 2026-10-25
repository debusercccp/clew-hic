COMPOSE := docker compose
LAB_PORT ?= $(shell grep -s '^LAB_PORT=' .env | cut -d= -f2)
LAB_PORT := $(or $(LAB_PORT),8888)
HIGLASS_PORT ?= $(shell grep -s '^HIGLASS_PORT=' .env | cut -d= -f2)
HIGLASS_PORT := $(or $(HIGLASS_PORT),8989)
NB ?= notebooks/01_explore.ipynb
LAB_PATH ?= lab

# Apre un URL nel browser se c'è un display grafico, altrimenti lo lascia stampato.
OPEN_URL = if [ -n "$$DISPLAY$$WAYLAND_DISPLAY" ] && command -v xdg-open >/dev/null; then \
	  xdg-open "$$url" >/dev/null 2>&1 & else echo "(nessun browser grafico: apri l'URL a mano)"; fi

.DEFAULT_GOAL := help
.PHONY: help env build up url notebook logs stop down restart test shell genome mcool higlass higlass-ingest clean distclean

help: ## Mostra i comandi disponibili
	@grep -E '^[a-z-]+:.*## ' $(MAKEFILE_LIST) | awk -F':.*## ' '{printf "  make %-15s %s\n", $$1, $$2}'

env: ## Crea .env da .env.example (UID/GID tuoi, password HiGlass casuale)
	@test -f .env || sed -e 's/^HOST_UID=.*/HOST_UID=$(shell id -u)/' \
	                     -e 's/^HOST_GID=.*/HOST_GID=$(shell id -g)/' .env.example > .env
	@grep -q '^HIGLASS_PASSWORD=.' .env || { grep -v '^HIGLASS_PASSWORD=' .env > .env.tmp; \
	  echo "HIGLASS_PASSWORD=$$(head -c 24 /dev/urandom | base64 | tr -dc 'A-Za-z0-9')" >> .env.tmp; \
	  mv .env.tmp .env; chmod 600 .env; }
	@mkdir -p results

build: env ## Costruisce l'immagine clew-env
	$(COMPOSE) build

up: env ## Avvia JupyterLab in background e stampa l'URL
	$(COMPOSE) up -d lab
	@$(MAKE) --no-print-directory url

url: ## Stampa l'URL di JupyterLab con il token
	@for i in $$(seq 30); do \
	  url=$$($(COMPOSE) exec -T lab jupyter server list 2>/dev/null | grep -o 'http://[^ ]*token=[a-z0-9]*' | head -1); \
	  if [ -n "$$url" ]; then echo "JupyterLab: $$(echo $$url | sed -E 's#^http://[^/]+/\??#http://127.0.0.1:$(LAB_PORT)/$(LAB_PATH)?#')"; exit 0; fi; \
	  sleep 1; \
	done; echo "JupyterLab non risponde: guarda 'make logs'" >&2; exit 1

notebook: env ## Avvia JupyterLab e apre il notebook nel browser (NB=notebooks/...)
	@test -f "$(NB)" || { echo "notebook non trovato: $(NB)" >&2; exit 1; }
	@$(COMPOSE) up -d lab
	@url=$$($(MAKE) --no-print-directory -s url LAB_PATH=lab/tree/$(NB)) || exit 1; \
	  url=$${url#JupyterLab: }; echo "Notebook: $$url"; $(OPEN_URL)

higlass: env ## Avvia il server HiGlass su localhost e crea l'utente per l'API
	@$(COMPOSE) --profile viz up -d higlass
	@for i in $$(seq 60); do curl -sf -o /dev/null http://127.0.0.1:$(HIGLASS_PORT)/api/v1/tilesets/ && break; sleep 2; done
	@$(COMPOSE) exec -T higlass python /home/higlass/projects/higlass-server/manage.py shell -c \
	  'import os; from django.contrib.auth.models import User; \
	   u, _ = User.objects.get_or_create(username=os.environ["HIGLASS_API_USER"]); \
	   u.set_password(os.environ["HIGLASS_API_PASSWORD"]); u.save()'
	@echo "HiGlass: http://127.0.0.1:$(HIGLASS_PORT)/app/  (registra i file con make higlass-ingest)"

higlass-ingest: higlass ## Registra .mcool in HiGlass e apre la vista (FILES="data/a.mcool ..." REGION=chr2:10Mb-20Mb ASM=hg38)
	@test -n "$(FILES)" || { echo 'specifica FILES, es. make higlass-ingest FILES=data/4DNFIXP4QG5B.mcool' >&2; exit 1; }
	@url=$$($(COMPOSE) run --rm -T worker python -m clew.ingest_higlass $(FILES) \
	    $(if $(REGION),--region $(REGION)) $(if $(ASM),--assembly $(ASM)) | tee /dev/stderr | grep -o 'http://[^ ]*config=[A-Za-z0-9_-]*'); \
	  test -n "$$url" || exit 1; $(OPEN_URL)

logs: ## Segue i log di JupyterLab (S=higlass per HiGlass)
	$(COMPOSE) --profile viz logs -f $(or $(S),lab)

stop: ## Ferma i container (JupyterLab e HiGlass) senza rimuoverli
	$(COMPOSE) --profile viz stop

restart: stop up ## Ferma e riavvia JupyterLab

down: ## Ferma e rimuove i container e la rete (immagini, dati e registrazioni HiGlass restano)
	$(COMPOSE) --profile batch --profile viz down --remove-orphans

test: env ## Esegue i test in un container worker usa-e-getta
	$(COMPOSE) run --rm worker

shell: env ## Shell nel container lab (lo avvia se serve)
	$(COMPOSE) up -d lab
	$(COMPOSE) exec lab bash

genome: ## Scarica il genoma di riferimento in data/genomes/ (ASM=hg38 | mm10)
	scripts/fetch_genome.sh $(or $(ASM),hg38)

mcool: ## Scarica un .mcool dal 4DN Data Portal in data/ (ACC=4DNFIXVAKX9Q)
	scripts/fetch_4dn.sh $(or $(ACC),$(error specifica ACC, es. make mcool ACC=4DNFIXVAKX9Q))

clean: down ## Come down, più l'immagine clew-env, la cache di build e il database HiGlass
	-docker volume rm clew_higlass-db
	-docker image rm clew-env:latest
	-docker builder prune -f --filter label=com.docker.compose.project=clew

distclean: clean ## Come clean, più cache Python e risultati in results/ (data/ non si tocca)
	find . -path ./data -prune -o \( -name __pycache__ -o -name .pytest_cache -o -name .ipynb_checkpoints \) -type d -print -exec rm -rf {} +
	rm -rf results
