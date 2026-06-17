## ---------------------------
##
## Script name: primary_B_bulk_rnaseq_processing.sh
##
## Original file: 260430.bulk.B.RNA.seq.processing.log.sh
##
## Purpose of script: Processing bulk B cell RNA-seq
##
##
## Date Created: 2026-04-30
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------




cd /path/to/output_directory

cat << 'EOF' > urls.txt
https://dl.dnanex.us/F/D/jZV4FZ4FBFj9fyPFJgQxf6JpzPkBQ2q7B49XFfy2/20260427_T4_72h_YT14040_S225_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/6pf5y95jQXFPJ7qGjq039Z571Bbxqb2BZ3fY1z4J/20260427_T7_24h_YT14040_S153_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/yqj2G8Vp3bjxFFKqV2jVzQbQ751ZG8ZFzBzZQv0x/20260427_T4_0h_YT14040_S118_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/Y0KVJxQb0J4Q6F7XgZPfb0Kqf2qXjj7Bx7V0K4Zx/20260427_C8_4h_YT14040_S130_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/7JgQ9V2ypB0JY5Yf5Z3jfqy002zQB5YZQXg8BVF5/20260427_C1_0h_YT14040_S107_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/0YYjYVF6f5VZBvbZ0KZJ8kPgYgqYQ45VqVY4kf4p/20260427_C5_4h_YT14040_S127_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/6776KZQ4v6QVj3GX8901701V12z54Xyx32fK1p20/20260427_C4_72h_YT14040_S158_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/QY05p5BZG0jq8z7jZG2XG0yK0PzgK8qz3xB4kPky/20260427_T2_72h_YT14040_S223_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/J8YXfb7x3XG62bpVf0xZ4vBbfzqkf1jJGYj9ZGV1/20260427_T4_72h_YT14040_S225_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/05x00Kz5yb206vjZygfy2Pj203Gv3Z4b68PXv7qX/20260427_T1_24h_YT14040_S147_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/kKyVPXp1kvQx40xjYPxx0FjK3FP9gPkP6yz5kB1v/20260427_T8_24h_YT14040_S154_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/7kbvZzZ1G0BB1p2PG4qg5y8pP6z9X91f152zpPJG/20260427_C5_72h_YT14040_S218_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/kzbVYvz1PGQBfpPQyQKq8P0j3vJV79Kv270y1PBK/20260427_T8_72h_YT14040_S229_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Qk5YJVYQjYx3PZzbFYf0XVk7Qp2Kxj482kfbjXby/20260427_C8_24h_YT14040_S146_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/33q4jQvZy16KPZ3GF4Px7YZ4GbkKx1xqpkFP5FBj/20260427_C6_24h_YT14040_S144_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/0Vgy2BQkg66ZXQK28fvFbQ5zp5jK6Y2Pb5bxbZVZ/20260427_T1_72h_YT14040_S222_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/y7qY5QB33ZKGFF9PzBy3j02Kz611gXv1pyJyg6Jb/20260427_T5_0h_YT14040_S119_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/2QQvF9ppgzBvzj8yzz8B9jj7BGQVPJKBZxq3V4bv/20260427_T7_72h_YT14040_S228_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/bYk3X97fKy26PZJvJp1FkBFqV1K1ZXz58xZPYFkZ/20260427_T8_72h_YT14040_S229_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/4VKPB3k3X8253bF5v3xzPjp8K183by9y71XQ9V9G/20260427_T5_72h_YT14040_S226_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/V7G4X2gk1yz39QPxjFKbk97K30ZQ92v6g877qf0f/20260427_C4_4h_YT14040_S126_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/bp654qYYgBK6PvV9z56ZGvJq77Kxg4XY0bzvKjGg/20260427_T7_0h_YT14040_S121_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/V1366k1yX5p42GQZ4YVZ81yy97JYjBYQxB8PK8jJ/20260427_C3_4h_YT14040_S125_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/qb8bb8xz2KYkzJvX63Z70vzzFFG8xY29B7Q15VqP/20260427_C2_0h_YT14040_S108_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/KPxv3YVqg4zVZpFzYqVx62Yp32GgfB20ZJ6b2vGZ/20260427_T1_4h_YT14040_S131_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/Y415g63pvGqxfbJbJqPj413zXxFzfj779Y3kgVXX/20260427_T6_0h_YT14040_S120_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/Gb9bzk2gPZvk85x60Y2xGY13b92Ypkz3FFk8bPVJ/20260427_T6_72h_YT14040_S227_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/vxQXZQkkXBBBYxP8yQxzGfQ11vbqb67P2X7zvjvp/20260427_T7_72h_YT14040_S228_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/8fBKpZ3zQf2bkKg2GPQY0K6vPgQ5ZK9Q39Xk014P/20260427_T2_24h_YT14040_S148_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/ZXJJJK0XVvVyZ29v4VzyPvpPJKKFyg04XByzpQYY/20260427_C2_72h_YT14040_S156_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/QjK86y8XQ5KZkY0Fx6YjZ858q09y6KyZ68gY3Z7J/20260427_C7_4h_YT14040_S129_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/2KZgxgy63Q3G690V55vqjj6jjFk7yFq9JZfGpJ1Z/20260427_C1_72h_YT14040_S155_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/gGG6kGfFYX7X49J1kBBxyb9y64kxfZJY1kzjqGQj/20260427_T3_4h_YT14040_S133_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/k02GjXk0qyKpQ16Z9YqGgxF46qZ5Fj69jVZ3YvKx/20260427_C4_0h_YT14040_S110_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/kz5Yff746234YVky08ZVBVX7pZ5PppvXyZ551YFF/20260427_C3_72h_YT14040_S157_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/B6pXzz26FPG06qg3g1Fxy3Vg711p0571qZygX6v2/20260427_T5_72h_YT14040_S226_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/G4Z4k4gKGPvXgQ522ZYjggg61zFBzQVYVpX1gBzy/20260427_C2_0h_YT14040_S108_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/PpQ4pbBvFb0GJyXyK9FP6VPBk5fkyf0J516zjf2V/20260427_T1_0h_YT14040_S115_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/KqQ5y89XzBJ1Vx95XYXfYB21fq30qVyjkqXf3pg5/20260427_C3_0h_YT14040_S109_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/JF15FPyzxX7G6j0bbVFvfbBG857bKjGBP7kzYFBY/20260427_C1_24h_YT14040_S139_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/39pvFqkyFQY3yYJZp95zF66pyq7FgbFfkfyZ2ZGF/20260427_C3_0h_YT14040_S109_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/8V206k59KY5Jy6jXP894qQXbB3Pbk22pxqB4v90X/20260427_C6_4h_YT14040_S128_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/zX8b3zpvFFZ81Py1YZfK1G31FYvGbQKvbKVPJ9vB/20260427_T3_24h_YT14040_S149_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/yjzjzKzv2Q68063kv22xf0zG9X9g1bKBFz4JQbBq/20260427_C3_4h_YT14040_S125_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/GGFXBB7byyy398qkGB636F26B29G19yX9F24jkZF/20260427_T3_72h_YT14040_S224_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/vYQVZ2YF50b70gk6yQf9z0j1Yyj4F8kg33KGqq77/20260427_C5_24h_YT14040_S143_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/4xVqvFyP5QQF9Y9Fvg27ZGJF67KX9BkZgQkB6P2z/20260427_T4_24h_YT14040_S150_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/0X9yQX9fz1jPqzZFX07QGqQF3vfFkvXJQb52FzvF/20260427_C5_72h_YT14040_S218_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/bgBFf7Gbzq9KpX1Z4FqyFQBGg9vGq4ygY2F0374B/20260427_C4_0h_YT14040_S110_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/xPgzV6Gyk1xQ0J6BqKqk7bp6Fg2J9fpx61pg4pzK/20260427_C4_24h_YT14040_S142_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/7ggp0z7v2pG215GYkZV33b28b1Jy11zXPPgP0vF5/20260427_T6_4h_YT14040_S136_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/pxpqJjPZq36yj2g24y7q4kKFvBZf07qYx3q6V2bb/20260427_C3_24h_YT14040_S141_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/9p4Q8Z2y5X4BgBjpPPZ65ZXxv6QVXbK72PxZBvVy/20260427_C6_72h_YT14040_S219_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/3P0yxk590gV3jbj3x1fG4BX5B5j2ybFzvv2KpJjZ/20260427_T1_72h_YT14040_S222_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/kxgvY933y0yq29qkvvXvb3q2F7B33VpgPGXpK4jV/20260427_C7_72h_YT14040_S220_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/PjgB0v2KxG62XPGyK5xFyX44ZfFVbZ4YPXvbKzJ7/20260427_T1_0h_YT14040_S115_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/FKGyjQJ03XJK33vgXZz6Q7z6fJ45JJ3X1fqPzjXg/20260427_C7_0h_YT14040_S113_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/41927vVx1y7bBxqX31g7Xz348Qx922JbP4bqKQ85/20260427_T8_0h_YT14040_S122_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/j97ZZyQq2XV9G86FF75j66PYx5P327x6fPfFyfVb/20260427_T3_0h_YT14040_S117_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/77359qpJkFXJy2y5xxYfkP5b24FyQQzx0vPJKKf8/20260427_T1_24h_YT14040_S147_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/V1x0PP5BfQKv31vZJKfG1j5g18YyKXYz8F86kbPB/20260427_C3_72h_YT14040_S157_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/qp449kYbGG2gFq2fP6Vvv54051jY6fV43zy5vFky/20260427_T5_0h_YT14040_S119_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/zjgZpk5KzGKzx2vYbQzp3x5J6Y5QQvPKfyJ9Q3V7/20260427_T8_4h_YT14040_S138_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/b83J45Y4ZX7kVKQP0bv4bxgg9kYVxGK80qJX01xF/20260427_T2_24h_YT14040_S148_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/0XV12K6vG44pqFgvgPYGZg9vK4QGj6Q87BJ0K7bX/20260427_T4_4h_YT14040_S134_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/xyYxz8x47f8fX1VQFjF95j86yg773fBQgY96KJ84/20260427_C6_24h_YT14040_S144_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Q83YV6vK5G05f0pVb2YxBPkq456gJfBZ49jQvp2J/20260427_T6_24h_YT14040_S152_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Kv4VY0JBjKX883Qk3bk55BGz3Fzp6PkG0B1JK82v/20260427_T5_24h_YT14040_S151_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/8kb0yGBXJ86X8X0QZgk5P2y0KYfQkXv8Kk2k8fpQ/20260427_C2_72h_YT14040_S156_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/YX8VB7zQPQ7V8KY2YJJyV35xyQzBVGgqv0q6pjk8/20260427_C6_4h_YT14040_S128_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/BJ9j9GQ59v2Qvpv6Z2Q892Bgz2Bv390zYKgGb1b6/20260427_C5_4h_YT14040_S127_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/vjJ9xp0yJqzKQqpB896B27qjjQKKpvkY26878QyK/20260427_T7_0h_YT14040_S121_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/5KB52Bq57Kv9ykPyxXXxB6b4J9BkfK92p0jgP6jq/20260427_C4_4h_YT14040_S126_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Z431VG7GgVJ6F3j95JYggP5FVgY21j69019J1p29/20260427_C1_0h_YT14040_S107_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/q4q5xbz0pY6j1gKGK8Y857Fpj3Q0PzGbFVgBfF0Y/20260427_C1_4h_YT14040_S123_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/98J14y7P7qKKgbYY79kgjy6P9qxGQy2KpGvfVbfQ/20260427_T6_4h_YT14040_S136_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/kz14pv03XXkbf82J8GP39vBv5pzjppPkVPz1Bv96/20260427_T3_0h_YT14040_S117_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/g3159pB57j6JZB0Z27ZB67663Q2FGJ5qKXYq4YbP/20260427_T4_24h_YT14040_S150_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/GqBGZpK2gQzvY8FvXFX5jy20201Z7vz2x6p8p3q4/20260427_T5_4h_YT14040_S135_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/k3bxpGgjjB4bZ3Z0gBXpbBf4G961k5091zpXyKq0/20260427_T2_4h_YT14040_S132_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/4PF6kx002XV0Bqg52yjYxYJX8Qx9FyqyJ3F6VVjz/20260427_C5_24h_YT14040_S143_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/6BB3bXYvP32yy0FzB3FJybJK9y2GzvbG8KxZ4V2V/20260427_C2_24h_YT14040_S140_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/BP713pyjgj701PFv991G9Yj8p7zpbKbY5KPPVZGJ/20260427_T5_4h_YT14040_S135_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/bxB13097Gjfk5KQk8qxzxx579VfbgPjzzQzG2qvj/20260427_T2_72h_YT14040_S223_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/x956Xf9ZPgQFK4qz9PyfY3Xj384k1VZpVg058kgb/20260427_T8_24h_YT14040_S154_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/F976VzzjfGZq03yjzjF19gYQBp8zQ0p9q0jXy6Jv/20260427_T2_0h_YT14040_S116_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Fq96V9KbxF31kfYGQ4gKBy5f9kB8By7VZqQ41YKB/20260427_C3_24h_YT14040_S141_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/1j6FVF2vQ8V00Jq62JG4g5B3V37QPK65P5Q6KQyG/20260427_T2_0h_YT14040_S116_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/49f6g9QV9kKbY98bFJ8jFyZYk98K1vK6Zf7q95jY/20260427_T3_4h_YT14040_S133_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/F7FKPK2BfQXjpfyqb9qjYVBFx3Y8g6b59XzfXpx6/20260427_C7_4h_YT14040_S129_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/84j2qb1Vz8G2FV6V591vZ3yPx0z1ZqkpYkgvjBqV/20260427_C6_0h_YT14040_S112_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/PGY5vj8KV1zYfg7z80XbyZzkXPXz0j9Zxz3yV9pY/20260427_C7_72h_YT14040_S220_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/89yX051kyk05qxZp7pYq9yP3QjK0FjvQ16x253yq/20260427_C1_4h_YT14040_S123_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/9PGVjvZ5V7X20q31pX0YFG40pp846q2b5QYykb1x/20260427_C2_4h_YT14040_S124_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/0B214fqyxpZF60K98x4zzqzBGPBFx3BbGZvQ7XZB/20260427_T7_24h_YT14040_S153_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/zBqgk344pyP11PJ9G8B6904pX1VyK94182Xf8F69/20260427_C7_24h_YT14040_S145_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/K4yBzbxv19vkJ33jYy2g2b74xxqxXF7VQ56v2y3g/20260427_C4_72h_YT14040_S158_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/FXbVvxJ5Bjz0ZzV09j4JkpbyyvzZ645299X45JKx/20260427_T4_4h_YT14040_S134_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/GfF0ff3jvPZk014Jv538QBKgzqzG59zZKPK2gvj9/20260427_T2_4h_YT14040_S132_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/P6308k1VB3J87VPxQzbVP84Z2F1jg7xFJbVY46fJ/20260427_T3_72h_YT14040_S224_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/2Kq6kZjyXXjk943p5B2P0P5gJGJpqPP0G12ZX17Y/20260427_C8_72h_YT14040_S221_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Vpkfgk78jP5kjzp6y9zvQY8z29g3bbxG3jGyQzv0/20260427_T7_4h_YT14040_S137_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/05614KK4B67pXJG1f3gGyzqxk8VbGXpbkQJbJqPk/20260427_T3_24h_YT14040_S149_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/JyvPJp39xzy6576Fy2kJg5J83fpVJyZJFFYgJ84Z/20260427_C4_24h_YT14040_S142_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/y4qp562q3X7ZJ42K7x8xJ2J577VXkv7bbyykp0y8/20260427_C1_72h_YT14040_S155_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/8VvXP95j8y78fKybyGgZby10xzxj448g66FK6PKy/20260427_T8_0h_YT14040_S122_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/7Gx3b7j0X21871BzKfJ1GjfgP60b87bj9y4qJbvk/20260427_C8_0h_YT14040_S114_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/PYbZ0Y791yYQFPBGBvYYXf9jFbpVxyYy2PqkPB5Q/20260427_T1_4h_YT14040_S131_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/ZFZy6b642BFXfk6z0y1z81X00BXX3X7Fby3bYY69/20260427_T8_4h_YT14040_S138_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/3qg8z8PzX66y0zkQV8Pyq31xzbF6pGX74yzxf93G/20260427_C8_24h_YT14040_S146_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/9qJgB3BJZgPxj6ZJFk7kGFf5GKKQq4Zj53K3J384/20260427_C5_0h_YT14040_S111_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/fG41p4VfP23Yz55ZG474Xj7gKqKJ14119zy84gK3/20260427_C6_0h_YT14040_S112_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/bpG2p1FvZ45889q4jJfvzxxzPZG8P5y448V449v8/20260427_C6_72h_YT14040_S219_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/B3fxy2x2P91qyqK2Z3bqgfXpVFY43Vy7VFy5831x/20260427_T7_4h_YT14040_S137_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/P8GZ7vjpFzQq8qGgpJ9zz5Fv26Gpb1k2QvJx4pjG/20260427_T5_24h_YT14040_S151_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/ykq8GGP7X4b1376zqqPVv2QbJXPv2gGyvxygY6qF/20260427_T4_0h_YT14040_S118_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/P371gxq165jQQ0b5868gb2B07k2gkgGfP8By10FX/20260427_C8_4h_YT14040_S130_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/ZX3xP7GJzP7Xpbq7bfGGZQPq8Bb4zzbF7BfY1Vjg/20260427_C5_0h_YT14040_S111_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Xb33zkXYkf4vkbzyXyvx69QyzJjK6x8ZZBXpKGqZ/20260427_C2_24h_YT14040_S140_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/fF6zy7vzzJqb29x6KGyVg95Z3q6yPK98kYFPqpgg/20260427_T6_0h_YT14040_S120_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/Fy0QYBXY4GkX57vVVjgV9GX529502zJj18vZ697V/20260427_C2_4h_YT14040_S124_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/PQJgZVfz976x4v4zXKKB9JV05Pg9bxjQX0vbVybj/20260427_T6_24h_YT14040_S152_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/gv3v5Qv8XpfGgjfFbxp07k9BZVyG5XJj42bF9bB1/20260427_C7_0h_YT14040_S113_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/21zPBKP44GV0v6pjQygk1jqBQQfPz6zg4fJ3zFKj/20260427_T6_72h_YT14040_S227_L008_R1_001.fastq.gz
https://dl.dnanex.us/F/D/z33QKxjp64BVfjq5yK40X576p47vj8kb2ZBP0VJP/20260427_C7_24h_YT14040_S145_L006_R2_001.fastq.gz
https://dl.dnanex.us/F/D/JqjxG5f3JY8JBY566J9bGg9b3f4Z1GYvXPJp4GZg/20260427_C8_72h_YT14040_S221_L008_R2_001.fastq.gz
https://dl.dnanex.us/F/D/ZYV5b5Q1K6Z8B489gYzgJB1qGGG1jP3qzGjjB8z8/20260427_C1_24h_YT14040_S139_L006_R1_001.fastq.gz
https://dl.dnanex.us/F/D/gbP2J7p2zGJg186xKxzVFF8V7jyFg05P5k95p0vb/20260427_C8_0h_YT14040_S114_L006_R2_001.fastq.gz
EOF


mkdir -p /path/to/output_directory
cd /path/to/output_directory



<<dl.code.sh
#!/bin/bash
#SBATCH --job-name=DL        # job name
#SBATCH --output=DL.out       # std out name
#SBATCH --error=DL.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=8         # num cpu cores / task
#SBATCH --mem=8G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short

THREADS=8
cd /path/to/output_directory
cat /path/to/urls.txt | xargs -P $THREADS -n 1 wget -c -t 5 --waitretry=5

if [ $? -eq 0 ]; then
    echo "All downloads completed successfully."
else
    echo "Some downloads might have failed. Please check the logs or re-run the script."
fi
dl.code.sh




#1. make layout file for analysis
cd /path/to/output_directory

for genotype in C T; do
for ID in `seq 1 8`; do
for time_p in 0h 4h 24h 72h; do
label="${genotype}${ID}_${time_p}"
R1_file=`ls *_${label}_*R1*.fastq.gz`
R2_file=`ls *_${label}_*R2*.fastq.gz`
echo -e "${label}\t/path/to/input_data/${R1_file}\t/path/to/input_data/${R2_file}"
done
done
done | sed '1i name\tR1_file\tR2_file' > layout.file.tsv



#2. run STAR mapping
mkdir -p /path/to/output_directory
cd /path/to/output_directory

awk 'NR>1 {print $0}' /path/to/layout.file.tsv | while read line; do
name=`echo -e ${line} | awk '{print $1}'`
R1_file=`echo -e ${line} | awk '{print $2}'`
R2_file=`echo -e ${line} | awk '{print $3}'`

if [ ! -d /path/to/output_directory ]; then
sbatch --nodelist=dn[028-029,031,034,038,040-045,052,062-063] /path/to/bulk_rnaseq_process_like_craft.sh \
/path/to/output_directory \
${R1_file} ${R2_file} ${name}
fi
done 






# run onslurm
cd /path/to/output_directory


<<code.sh
#!/bin/bash
#SBATCH --job-name=summary         # job name
#SBATCH --output=summary.out       # std out name
#SBATCH --error=summary.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=1         # num cpu cores / task
#SBATCH --mem=32G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short

export OMP_NUM_THREADS=1 
eval "$(conda shell.bash hook)"


#3. merge into text file
mkdir -p /path/to/output_directory
cd /path/to/output_directory


conda activate R
Rscript -e '
library(tidyverse)
gene_file="/path/to/geneInfo.tab"
layout_file="/path/to/layout.file.tsv"

layout_df<-read_tsv(layout_file)
for(i in 1:nrow(layout_df)){
	sample_name<-layout_df$name[i]
	file_name<-str_c("/path/to/output_directory",sample_name,"/STARoutReadsPerGene.out.tab")
	tmp<-read_tsv(file_name,col_names=F)%>%dplyr::select(X1,X2)
	colnames(tmp)<-c("gene_id",sample_name)
	if(i==1){
		result<-tmp
	}else{
		result<-full_join(result,tmp,by="gene_id")
	}
}

gene_df<-read_tsv(gene_file,skip=1,col_names=F)
colnames(gene_df)<-c("gene_id","Gene","category")

gene_result<-inner_join(gene_df,result,by="gene_id")

write_tsv(gene_result,"Bulk.CD40L.Daudi.RNA.seq.count.annotation.tsv.gz")
write_tsv(result,"Bulk.CD40L.Daudi.RNA.seq.count.tsv.gz")
'



Rscript -e '
library(tidyverse)
gene_file="/path/to/geneInfo.tab"
layout_file="/path/to/layout.file.tsv"

layout_df<-read_tsv(layout_file)
layout_df<-filter(layout_df,!str_detect(name,"C6"))
for(i in 1:nrow(layout_df)){
	sample_name<-layout_df$name[i]
	file_name<-str_c("/path/to/output_directory",sample_name,"/STARoutReadsPerGene.out.tab")
	tmp<-read_tsv(file_name,col_names=F)%>%dplyr::select(X1,X2)
	colnames(tmp)<-c("gene_id",sample_name)
	if(i==1){
		result<-tmp
	}else{
		result<-full_join(result,tmp,by="gene_id")
	}
}

gene_df<-read_tsv(gene_file,skip=1,col_names=F)
colnames(gene_df)<-c("gene_id","Gene","category")

gene_result<-inner_join(gene_df,result,by="gene_id")

write_tsv(gene_result,"Bulk.B.cell.remove.C6.Daudi.RNA.seq.count.annotation.tsv.gz")
write_tsv(result,"Bulk.B.cell.remove.C6.CD40L.Daudi.RNA.seq.count.tsv.gz")
'





Rscript -e '
library(tidyverse)
gene_file="/path/to/geneInfo.tab"
layout_file="/path/to/layout.file.tsv"

layout_df<-read_tsv(layout_file)
for(i in 1:nrow(layout_df)){
	sample_name<-layout_df$name[i]
	file_name<-str_c("/path/to/output_directory",sample_name,"/STARoutReadsPerGene.out.tab")
	tmp<-read_tsv(file_name,col_names=F)%>%dplyr::select(X1,X2)
	colnames(tmp)<-c("gene_id",sample_name)
	if(i==1){
		result<-tmp
	}else{
		result<-full_join(result,tmp,by="gene_id")
	}
}

gene_df<-read_tsv(gene_file,skip=1,col_names=F)
colnames(gene_df)<-c("gene_id","Gene","category")
colnames(result)<-gsub("C6_","T9_",colnames(result))
gene_result<-inner_join(gene_df,result,by="gene_id")

write_tsv(gene_result,"Bulk.B.cell.C6.to.T9.Daudi.RNA.seq.count.annotation.tsv.gz")
write_tsv(result,"Bulk.B.cell.C6.to.T9.CD40L.Daudi.RNA.seq.count.tsv.gz")
'





#4. Summarize fastQC
conda activate multiqc-env  
cd /path/to/output_directory

BASE_DIR="/path/to/output_directory"
OUT_DIR="${BASE_DIR}/summary/raw.multiQC_results"

for i in 1 2; do
  TMP_DIR="/path/to/output_directory"
  mkdir -p ${TMP_DIR}

  echo "Finding FastQC files for R${i}..."
  find ${BASE_DIR} -type f -name "20260427_*_R${i}_*_fastqc.zip" -exec cp {} ${TMP_DIR}/ \;

  echo "Running MultiQC for R${i}..."
  multiqc ${TMP_DIR} \
    -n "multiqc_report_pilot_R${i}" \
    --outdir ${OUT_DIR} \
    --title "Pilot_Sequencing_R${i}"

  rm -r ${TMP_DIR}
done


OUT_DIR="${BASE_DIR}/summary/trimmed.multiQC_results"
for i in 1 2; do
  TMP_DIR="/path/to/output_directory"
  mkdir -p ${TMP_DIR}

  echo "Finding FastQC files for R${i}..."
  find ${BASE_DIR} -type f -name "trimmed*.R${i}_*fastqc.zip" -exec cp {} ${TMP_DIR}/ \;

  echo "Running MultiQC for R${i}..."
  multiqc ${TMP_DIR} \
    -n "multiqc_report_pilot_R${i}" \
    --outdir ${OUT_DIR} \
    --title "Pilot_Sequencing_R${i}"

  rm -r ${TMP_DIR}
done

code.sh




#for sanity check, call genotypes from bam files
mkdir -p /path/to/output_directory
cd /path/to/output_directory

awk 'NR>1 {print $0}' /path/to/layout.file.tsv | while read line; do
name=`echo -e ${line} | awk '{print $1}'`
if [ ! -d /path/to/output_directory ]; then
sbatch --exclude=dn054 \
/path/to/call_genotype_from_bam.sh \
/path/to/output_directory \
/path/to/STARoutAligned.sortedByCoord.out.bam \
${name}
fi
done
