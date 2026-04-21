# Development image for kronos-semi.
#
# Based on the official dolfinx image, which provides FEniCSx (dolfinx, ufl,
# basix, ffcx), PETSc, MPI, and a compatible Python. We layer the project's
# Python dev dependencies and install the package in editable mode on top.
#
# Build:   docker build -t kronos-semi:dev .
# Run:     docker compose up -d dev   (preferred, see docker-compose.yml)

FROM ghcr.io/fenics/dolfinx/dolfinx:stable

ARG USERNAME=dev
ARG USER_UID=1000
ARG USER_GID=1000

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# System tools useful inside a dev container. The base image already has
# build-essential, git, and Python; we add interactive conveniences.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        sudo \
        vim \
        less \
        curl \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Create a non-root user that matches the host UID/GID so bind-mounted files
# stay owned by the host user. If a user or group already exists at the
# requested UID/GID (the dolfinx base image ships an `ubuntu` user at 1000),
# rename it to ${USERNAME} rather than failing.
RUN set -eux; \
    if getent group ${USER_GID} >/dev/null; then \
        existing_group="$(getent group ${USER_GID} | cut -d: -f1)"; \
        if [ "${existing_group}" != "${USERNAME}" ]; then \
            groupmod -n ${USERNAME} "${existing_group}"; \
        fi; \
    else \
        groupadd --gid ${USER_GID} ${USERNAME}; \
    fi; \
    if getent passwd ${USER_UID} >/dev/null; then \
        existing_user="$(getent passwd ${USER_UID} | cut -d: -f1)"; \
        if [ "${existing_user}" != "${USERNAME}" ]; then \
            usermod -l ${USERNAME} -d /home/${USERNAME} -m "${existing_user}"; \
            usermod -g ${USER_GID} ${USERNAME}; \
        fi; \
    else \
        useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME}; \
    fi; \
    echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME}; \
    chmod 0440 /etc/sudoers.d/${USERNAME}

WORKDIR /workspaces/kronos-semi

# Install Python dev dependencies first so the layer caches across code edits.
# Using --system-site-packages semantics: dolfinx is already importable from
# the base image's Python; we install into the same interpreter.
COPY pyproject.toml README.md ./
RUN pip install --break-system-packages \
        "numpy>=1.24" \
        "matplotlib>=3.7" \
        "jsonschema>=4.0" \
        "pytest>=7.0" \
        "pytest-cov>=4.0" \
        "nbformat>=5.0" \
        "ruff>=0.1" \
        "jupyterlab>=4.0" \
        "ipykernel>=6.0"

# Copy the rest of the source and install the package editable. The bind
# mount in docker-compose will overlay /workspaces/kronos-semi at runtime,
# but installing here ensures `import semi` works if the image is run without
# a bind mount (e.g., CI).
COPY . .
RUN pip install --break-system-packages -e ".[dev]"

RUN chown -R ${USER_UID}:${USER_GID} /workspaces

USER ${USERNAME}

CMD ["bash"]
