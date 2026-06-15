# antoinejhaddad.com

A single-page personal / bioinformatics site. Plain HTML + CSS, no build step,
hosted for free on GitHub Pages with the custom domain **antoinejhaddad.com**.

```
index.html      the page
styles.css      all styling (palette pulled from a eucalyptus trunk photo)
favicon.svg     four-colour palette mark
CNAME           custom domain for GitHub Pages
.nojekyll       tells GitHub Pages to serve files as-is
```

---

## 1. Fill in your details (2 minutes)

Open `index.html` and replace these placeholders (search for `EDIT`):

- **LinkedIn** — `https://www.linkedin.com/in/your-profile` (appears twice)
- **Instagram** — `https://instagram.com/your-handle` (appears twice)
- **GitHub** — `https://github.com/your-username` (appears twice)
- **Education / Focus / Projects** lines in the Bioinformatics section
  (the `[University]`, `[...]` bits)

Your email `connect@antoinejhaddad.com` is already wired up.

## 2. Put it on GitHub

You get GitHub Pages free as a student (and the
[Student Developer Pack](https://education.github.com/pack) on top).

```bash
# from inside this folder
git init
git add .
git commit -m "Personal site"
git branch -M main

# create a new EMPTY repo on github.com named:  Website
# (no README / .gitignore / license), then:
git remote add origin https://github.com/arangutambo/Website.git
git push -u origin main
```

## 3. Turn on GitHub Pages

On GitHub: **Settings → Pages**

- **Source:** Deploy from a branch
- **Branch:** `main`  ·  **Folder:** `/ (root)`  →  Save

In a minute your site is live at `https://arangutambo.github.io/Website/`
(and at `https://antoinejhaddad.com` once DNS is set up in step 4).

## 4. Point your domain at it

The `CNAME` file already contains `antoinejhaddad.com`, so GitHub will pick it up.
At your domain registrar (where you bought antoinejhaddad.com), add these DNS records:

**Apex domain (antoinejhaddad.com) — four A records:**

```
A   @   185.199.108.153
A   @   185.199.109.153
A   @   185.199.110.153
A   @   185.199.111.153
```

(Optional, recommended — IPv6 AAAA records:)

```
AAAA  @  2606:50c0:8000::153
AAAA  @  2606:50c0:8001::153
AAAA  @  2606:50c0:8002::153
AAAA  @  2606:50c0:8003::153
```

**www subdomain — one CNAME record:**

```
CNAME   www   arangutambo.github.io.
```

Then back in **Settings → Pages**, set **Custom domain** to `antoinejhaddad.com`
and tick **Enforce HTTPS** once the certificate is issued (can take up to ~24h
after DNS propagates).

---

## Editing later

It's just one HTML file and one CSS file — open them in any editor (Obsidian
works fine), commit, and `git push`. Pages redeploys automatically.

The colour palette lives at the top of `styles.css` as CSS variables
(`--brown`, `--rust`, `--green`, `--purple`, `--paper`) if you ever want to
tweak the tones.
