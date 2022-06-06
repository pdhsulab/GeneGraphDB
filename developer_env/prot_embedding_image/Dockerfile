FROM gcr.io/durrant/genegraphdb_base:latest

# install pytorch for expected GPU version
RUN pip install --upgrade pip && \
    pip install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu113 && \
    rm -r ~/.cache/pip
