# How to Upload ChemMind to GitHub

Since this is a local project, you need to create a repository on GitHub and push your code to it.

## Step 1: Create a New Repository
1.  Go to **[https://github.com/new](https://github.com/new)**.
2.  Enter a repository name (e.g., `chemmind-pro`).
3.  **Do not** initialize with README, .gitignore, or License (we already have them).
4.  Click **Create repository**.

## Step 2: Push Your Code
Copy the commands from the "…or push an existing repository from the command line" section, or use the ones below.

Run these commands in your terminal (inside the `chemmind` folder):

```bash
# 1. Link your local repo to the remote GitHub repo
git remote add origin https://github.com/YOUR_USERNAME/chemmind-pro.git

# 2. Rename the branch to main (if not already)
git branch -M main

# 3. Push your code
git push -u origin main
```

> **Note**: You will need to replace `YOUR_USERNAME` and `chemmind-pro` with your actual GitHub username and repository name.
