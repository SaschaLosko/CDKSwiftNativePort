# Publishing

## 1) Create remote repository

Example (GitHub CLI):

```bash
gh repo create <your-org-or-user>/CDKPortExperimental --public --source . --remote origin --push
```

Or add remote manually:

```bash
git remote add origin git@github.com:<your-org-or-user>/CDKPortExperimental.git
```

## 2) Push main branch

```bash
git push -u origin main
```

## 3) Push tags

```bash
git push origin --tags
```

## 4) Consume from app

In Xcode package dependencies, point to:

- `https://github.com/<your-org-or-user>/CDKPortExperimental.git`

and choose the published tag (for example `0.1.0`).
