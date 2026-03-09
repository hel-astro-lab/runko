"""Cross-platform CPU RAM usage helper."""

import sys
import resource


def get_rss_kB() -> int:
    """Return the Resident Set Size of this process in kB.

    - Linux: reads VmRSS from /proc/self/status (current RSS).
    - macOS: resource.getrusage ru_maxrss in bytes (peak RSS).
    """

    if sys.platform == "linux":
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])  # already in kB
        # fallback if VmRSS line absent
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    # macOS (darwin) returns bytes; convert to kB
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
