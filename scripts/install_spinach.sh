#!/usr/bin/env bash
# Download the pinned Spinach release into this repository's ignored lib/ folder.
set -euo pipefail

repo_root="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
version="${SPINACH_VERSION:-2.10.1}"
target="${repo_root}/lib/Spinach-${version}"
archive_url="https://github.com/IlyaKuprov/Spinach/archive/refs/tags/${version}.zip"

if [ -d "${target}/kernel" ]; then
    echo "Spinach ${version} is already installed at ${target}"
    exit 0
fi

if [ -e "${target}" ]; then
    echo "Found an incomplete installation at ${target}." >&2
    echo "Remove it manually, then run this installer again." >&2
    exit 1
fi

command -v curl >/dev/null 2>&1 || {
    echo "curl is required to download Spinach." >&2
    exit 1
}
command -v unzip >/dev/null 2>&1 || {
    echo "unzip is required to install Spinach." >&2
    exit 1
}

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/spinach.XXXXXX")"
trap 'rm -rf "${tmp_dir}"' EXIT

echo "Downloading Spinach ${version}..."
curl --fail --location --retry 3 --output "${tmp_dir}/spinach.zip" "${archive_url}"
unzip -q "${tmp_dir}/spinach.zip" -d "${tmp_dir}"

source_dir="${tmp_dir}/Spinach-${version}"
if [ ! -d "${source_dir}/kernel" ]; then
    echo "The downloaded archive did not contain the expected Spinach kernel." >&2
    exit 1
fi

mkdir -p "${repo_root}/lib"
mv "${source_dir}" "${target}"

for required_dir in kernel etc experiments interfaces; do
    if [ ! -d "${target}/${required_dir}" ]; then
        echo "Spinach installation is missing ${required_dir}/" >&2
        exit 1
    fi
done

echo "Spinach ${version} installed at ${target}"
echo "Next: source ${repo_root}/setup_paths.sh"
