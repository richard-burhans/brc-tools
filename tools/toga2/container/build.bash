#!/usr/bin/env bash
# Build and optionally publish the TOGA2 images, one per upstream tag.
#
# ⛔ NO CREDENTIAL IS EVER WRITTEN IN THIS FILE OR PASSED ON A COMMAND LINE. `docker login
# --password <secret>` puts the secret in the process table, in shell history, and in any captured
# build log -- Docker itself warns about it on every invocation. This reads from the environment
# and pipes through stdin instead.
#
#   QUAY_ROBOT_USER=galaxy+toga2 QUAY_ROBOT_PASSWORD=... ./build.bash --push
#   ./build.bash --tags v2.0.9             # build only, no registry involved
set -o errexit -o nounset -o pipefail

REGISTRY="${REGISTRY:-quay.io/galaxy/toga2}"
# ⚠ DERIVE THE LOGIN HOST FROM $REGISTRY. These were previously independent: the host was
# hardcoded `quay.io` while the push used $REGISTRY, so `REGISTRY=ghcr.io/org/toga2 ./build.bash
# --push` authenticated to quay.io and then pushed to ghcr.io unauthenticated.
REGISTRY_HOST="${REGISTRY%%/*}"
REPO="https://github.com/hillerlab/TOGA2.git"
PUSH=0
TAGS=""

while [ $# -gt 0 ]; do
    case "$1" in
        --push)  PUSH=1 ;;
        --tags)  shift; [ $# -gt 0 ] || { echo "usage: $0 [--push] [--tags 'v2.0.8 v2.0.9']" >&2; exit 2; }
                 TAGS="$1" ;;
        *) echo "usage: $0 [--push] [--tags 'v2.0.8 v2.0.9']" >&2; exit 2 ;;
    esac
    shift
done

# ⚠ TAGS ARE DISCOVERED, NOT MAINTAINED IN A LIST HERE. A hand-kept list goes stale silently and
# the first sign is an image nobody built. Restrict with --tags when iterating.
#
# ⛔ BUT DISCOVERY WITHOUT A SHAPE FILTER IS THE OTHER FAILURE. Unfiltered, the day upstream adds
# a `latest` or `stable` ref this script silently redefines $REGISTRY:latest as whatever commit
# that ref points at. The grep accepts only vN.N.N with an optional letter suffix (v2.0.7a), which
# excludes moving refs and release candidates. `-f3-` rather than `-f3` because a tag containing a
# slash (release/2.2) would otherwise truncate to `release` and collide with its siblings.
if [ -z "$TAGS" ]; then
    TAGS=$(git ls-remote --tags --refs "$REPO" \
             | cut -d/ -f3- \
             | grep -E '^v[0-9]+\.[0-9]+\.[0-9]+[a-z]?$' \
             | sort -V | tr '\n' ' ')
fi

# ⚠ NO PATHNAME EXPANSION ON THE TAG LIST. Unquoted `$TAGS` in a for-loop is subject to globbing,
# not just word splitting, so a tag containing * ? or [ would silently become whatever matches in
# the current directory. Split once, explicitly, into an array.
set -o noglob
read -r -a TAG_LIST <<< "$TAGS"
set +o noglob
[ "${#TAG_LIST[@]}" -gt 0 ] || { echo "no tags matched" >&2; exit 1; }

if [ "$PUSH" -eq 1 ]; then
    # ⛔ TRACE OFF ACROSS THE WHOLE CREDENTIAL REGION. --password-stdin keeps the secret out of the
    # process table, but it does nothing about `set -x`: under xtrace both the :? guard and the
    # printf expand the token into stderr, which is exactly what a CI runner captures into a build
    # log. And xtrace can be switched on from OUTSIDE via an exported SHELLOPTS, so it is not
    # enough to simply not type -x. The redirect on `set +x` suppresses the trace of that line.
    _xtrace=0; case "$-" in *x*) _xtrace=1 ;; esac
    { set +x; } 2>/dev/null

    # ⛔ A PRIVATE, DISPOSABLE DOCKER CONFIG -- NOT ~/.docker/config.json. Without this, `docker
    # login` writes the robot token to the user's real config as plain base64 (no credential
    # helper), where it persists for anything that can read $HOME. Cleanup used to rely on
    # `docker logout`, which is unverifiable: it prints "Removing login credentials", exits 0, and
    # on a read-only config dir leaves the file byte-identical with the live secret still in it.
    # `rm -rf` of a directory we created is enforceable, and it is also per-invocation, so this no
    # longer destroys a pre-existing session or races a concurrent build.
    #
    # ⚠ THE TRAP IS REGISTERED BEFORE THE LOGIN. Registered after, a signal in the window between
    # them leaves the credential behind with no cleanup at all.
    DOCKER_CONFIG="$(mktemp -d)"
    export DOCKER_CONFIG
    trap 'rm -rf "$DOCKER_CONFIG"' EXIT

    : "${QUAY_ROBOT_USER:?set QUAY_ROBOT_USER}"
    : "${QUAY_ROBOT_PASSWORD:?set QUAY_ROBOT_PASSWORD}"
    printf '%s' "$QUAY_ROBOT_PASSWORD" \
        | docker login --username "$QUAY_ROBOT_USER" --password-stdin "$REGISTRY_HOST"

    [ "$_xtrace" -eq 1 ] && set -x
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PUSHED=()

for tag in "${TAG_LIST[@]}"; do
    version="${tag#v}"
    echo "== $tag -> $REGISTRY:$version"

    # ⛔ THE BUILD IS THE TEST -- and it is only a real test because the Dockerfile ends with a
    # complete synthetic TOGA2 run, not because it runs check_dependencies.py. That checker
    # verifies the PRESENCE of three third-party programs and passes on an image with Nextflow,
    # toga2.py and CESAR2.0 deleted. Do not add `|| true` to this line: a container that builds
    # and then emits a differently shaped orthology_classification.tsv is worse than no container.
    # ⚠ CONTEXT IS THE TOOL DIRECTORY, NOT container/. The Dockerfile's smoke test COPYs
    # test-data/, and those fixtures are the SAME files Galaxy's <test> in toga2.xml uses. Keeping
    # one copy at the Galaxy-conventional location and reaching up for it beats duplicating them.
    docker buildx build --progress=plain \
        --file "$here/Dockerfile" \
        --build-arg "TOGA2_TAG=$tag" \
        --tag "$REGISTRY:$version" \
        "$here/.."

    if [ "$PUSH" -eq 1 ]; then
        # ⛔ TAKE THE DIGEST FROM THE PUSH, NOT FROM `docker inspect`. Two things were wrong with
        # inspecting: `{{index .RepoDigests 0}}` returns whatever is FIRST, which for an image
        # sharing content with another tagged image is a digest for a repository this script never
        # pushed to; and a freshly built, never-pushed image already HAS a RepoDigest, so the line
        # printed a plausible repo@sha256 even when the push had failed. `docker push` reports the
        # digest the registry actually stored, which is the only one worth recording.
        digest=$(docker push "$REGISTRY:$version" | sed -n 's/.*\(digest: \)\?\(sha256:[0-9a-f]\{64\}\).*/\2/p' | tail -1)
        if [ -z "$digest" ]; then
            echo "ERROR: push reported no digest for $REGISTRY:$version" >&2
            exit 1
        fi
        echo "   pushed: $REGISTRY@$digest"
        PUSHED+=("$REGISTRY:$version -> $REGISTRY@$digest")
    fi
done

# ⚠ SAY WHAT ACTUALLY SHIPPED. Without a completion marker a run that died at tag 2 of 4 and a run
# that finished all 4 are distinguishable only by counting lines against a tag list you have to
# re-derive. Pin the wrapper to a digest from this list, never to the tag.
echo
echo "built: ${TAG_LIST[*]}"
if [ "$PUSH" -eq 1 ]; then
    printf 'pushed %d image(s):\n' "${#PUSHED[@]}"
    printf '   %s\n' "${PUSHED[@]}"
fi
