# Quick Start Guide

## Installation & Setup

1. **Install dependencies:**
   ```bash
   cd frontend
   npm install
   ```

2. **Run the development server:**
   ```bash
   npm run dev
   ```

3. **Open your browser:**
   Navigate to [http://localhost:3000](http://localhost:3000)

## What You'll See

- **Left Side**: Logo with clickable instructions
- **Right Side**: Prompt area for experiment creation

## Optional: Display a Walkthrough Video

1. Place your MP4 file in `frontend/public/` and name it `simulation-demo.mp4`
2. (Optional) Add a poster image `video-placeholder.png` in the same folder
3. After you run a simulation, the video will appear on the left panel in place of the logo

## Making Changes

### To modify the layout:
- Edit `app/page.tsx` for the main layout structure
- Edit `components/LogoSection.tsx` for the left side
- Edit `components/PromptSection.tsx` for the right side

### To change colors:
- Edit CSS variables in `app/globals.css`

### To add components:
```bash
npx shadcn-ui@latest add [component-name]
```

## Next Steps

1. Connect the frontend to the backend API
2. Implement result display after simulation
3. Add file download functionality
4. Enhance the UI with more features

## Troubleshooting

- **Port already in use?** Change the port in `package.json` scripts or use `-p 3001`
- **Build errors?** Make sure all dependencies are installed: `npm install`
- **Type errors?** Run `npm run build` to check for TypeScript errors

