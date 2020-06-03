#!/usr/bin/env sh

PROJECT='hw2'
PORT=8888
# HASH=`python -c "import IPython; print(IPython.lib.passwd('${PROJECT}'))"` # IPython (dependancy of Jupyter) must be install globally

# Pull/Build Docker image 
# docker pull jupyter/minimal-notebook
docker build . -t "${PROJECT}"

# Sleep for server to start up then open browser to URL (sent to background)
#sleep 5 && python -m webbrowser -n "http://127.0.0.1:8888" &
sleep 5 && open "http://127.0.0.1:${PORT}" &

# Run Docker image
#docker run --rm -it --name "${PROJECT}" -p "${PORT}":"${PORT}" -e JUPYTER_ENABLE_LAB=yes -v "${PWD}":/home/jovyan/work "${PROJECT}" start-notebook.sh --NotebookApp.password="${HASH}" # https://jupyter-docker-stacks.readthedocs.io/en/latest/using/common.html#notebook-options
docker run --rm -it --name "${PROJECT}" -p "${PORT}":"${PORT}" -v "${PWD}":/home/jovyan/work "${PROJECT}" start.sh jupyter lab --LabApp.token='' # https://jupyter-docker-stacks.readthedocs.io/en/latest/using/recipes.html
