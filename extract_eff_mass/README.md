# Extract mass

Script to fetch data from the Supercomputing Wales server and
generate the effective mass solving the cosh equation. Provided a
channel name, it downloads the suitable data from the server and
process it for all temperatures available, all flavor structures and
all types of sources available. In order to run the script, use

```bash
bash fetch_data.sh name_channel
```
where `name_channel` can be `g5, vec, ax_plus, ax_minus, g0`. The
produced is stored in a folder and separated into flavors ->
temperature -> source.


