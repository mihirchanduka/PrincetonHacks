# Virtual Mouse Lab - Frontend

A modern web interface for the Virtual Mouse Lab platform, built with Next.js, TypeScript, and shadcn/ui.

## Features

- 🎨 Research-inspired gray and beige color scheme
- 📱 Responsive split-screen layout
- 🖱️ Interactive logo with instructions
- 💬 Prompt-based experiment creation
- 📥 Download simulation results as ZIP files
- ⚡ Fast development with Next.js

## Getting Started

### Prerequisites

- Node.js 18+ and npm/yarn/pnpm
- Backend API running on `http://localhost:8000` (see main project README)

### Installation

1. Navigate to the frontend directory:
```bash
cd frontend
```

2. Install dependencies:
```bash
npm install
# or
yarn install
# or
pnpm install
```

3. Configure API URL (optional):
   - The default API URL is `http://localhost:8000`
   - You can override this by creating a `.env.local` file with:
     ```
     NEXT_PUBLIC_API_URL=http://your-api-url:8000
     ```

4. Run the development server:
```bash
npm run dev
# or
yarn dev
# or
pnpm dev
```

5. Open [http://localhost:3000](http://localhost:3000) in your browser

## Using the Application

1. **Enter a prompt**: Describe your experiment in the text area on the right side
   - Example: "Create a C57BL/6J mouse with Trp53 knockout. Apply NewDrugX at 20 mg/kg for 30 days. Run blood work and open field test."

2. **Click "Run Simulation"**: The system will parse your prompt and run the simulation

3. **Download Results**: Once complete, click "Download Results (ZIP)" to get the simulation package

## Prompt Examples

The system can understand simple prompts and extract:
- **Genotype**: Mouse strain (C57BL/6J, BALB/c) and knockouts (Trp53)
- **Treatment**: Drug name and dosage (e.g., "NewDrugX at 20 mg/kg")
- **Tests**: Blood work, open field test, behavioral tests

### Example Prompts:

- "Create a C57BL/6J mouse with no modifications. Run blood work and open field test."
- "Generate a mouse with Trp53 knockout. Apply NewDrugX at 20 mg/kg for 30 days. Run all tests."
- "Create a mouse and simulate high-fat diet intervention. Analyze behavioral changes and biomarkers."

## Project Structure

```
frontend/
├── app/
│   ├── globals.css          # Global styles and color scheme
│   ├── layout.tsx           # Root layout
│   └── page.tsx             # Main page
├── components/
│   ├── ui/                  # shadcn/ui components
│   ├── LogoSection.tsx      # Left side logo and instructions
│   └── PromptSection.tsx    # Right side prompt area
├── lib/
│   └── utils.ts             # Utility functions
└── public/                  # Static assets
```

## Color Scheme

The application uses a research-inspired color palette:
- **Background**: Light brown/tan (`#D4C4B0`)
- **Foreground**: Dark brown text (`#3A2E1F`)
- **Primary**: Darker brown (`#5A4A3A`)
- **Borders**: Light blue (`#7A9CC6`)
- **Accent**: Medium brown (`#B8A088`)

## API Integration

The frontend communicates with the backend API at `http://localhost:8000`:

- **POST /api/simulate**: Run a simulation
- **GET /api/download/{simulation_id}**: Download simulation results

Make sure the backend is running before using the frontend.

## Customization

### Adding shadcn/ui Components

To add more shadcn/ui components:

```bash
npx shadcn-ui@latest add [component-name]
```

### Modifying Colors

Edit the CSS variables in `app/globals.css` to customize the color scheme.

### Changing API URL

Create or edit `.env.local`:

```
NEXT_PUBLIC_API_URL=http://your-api-url:8000
```

## Building for Production

```bash
npm run build
npm start
```

## Technologies

- **Next.js 14** - React framework
- **TypeScript** - Type safety
- **Tailwind CSS** - Styling
- **shadcn/ui** - UI components
- **Radix UI** - Accessible component primitives
- **Lucide React** - Icons

## Troubleshooting

### CORS Errors

If you see CORS errors, make sure:
1. The backend API is running on `http://localhost:8000`
2. CORS middleware is enabled in `main_api.py`
3. The frontend is running on `http://localhost:3000`

### API Connection Issues

- Check that the backend is running: `uvicorn main_api:app --reload --port 8000`
- Verify the API URL in `.env.local`
- Check browser console for error messages

### Build Errors

- Make sure all dependencies are installed: `npm install`
- Clear `.next` folder and rebuild: `rm -rf .next && npm run build`
- Check for TypeScript errors: `npm run build`
