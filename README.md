![](plots/solution.png)

### Notices

Git mirrors:

- [Codeberg](https://github.com/paveloom-university/Celestial-Mechanics-Laboratory-Workshop-S09-2021)
- [GitHub](https://gitlab.com/paveloom-g/university/s09-2021/celestial-mechanics-laboratory-workshop-s09-2021)
- [GitLab](https://codeberg.org/paveloom-university/Celestial-Mechanics-Laboratory-Workshop-S09-2021)

This project provides [Julia](https://julialang.org) scripts. Make sure to use the project files (`Project.toml`) when running them:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/script.jl
```

Alternatively, you can use the `julia.bash` script, which starts a [daemon](https://github.com/dmolina/DaemonMode.jl) and runs scripts through it:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
./julia.bash scripts/script.jl
```

To kill the daemon run

```bash
./julia.bash kill
```
