"use client"

import { useEffect, useState } from "react"
import Image from "next/image"
import { Button } from "@/components/ui/button"
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "@/components/ui/dialog"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Beaker, Mouse, Dna, BarChart3 } from "lucide-react"

export function LogoSection() {
  const [isDialogOpen, setIsDialogOpen] = useState(false)
  const [showVideo, setShowVideo] = useState(false)
  const [videoKey, setVideoKey] = useState(0)

  useEffect(() => {
    if (typeof window === "undefined") return

    const handleToggle = (event: Event) => {
      const customEvent = event as CustomEvent<{ show: boolean }>
      const shouldShow = Boolean(customEvent.detail?.show)
      setShowVideo(shouldShow)
      if (shouldShow) {
        setVideoKey((prev) => prev + 1)
      }
    }

    window.addEventListener("simulation-video", handleToggle)

    return () => {
      window.removeEventListener("simulation-video", handleToggle)
    }
  }, [])

  return (
    <div className="flex flex-col items-center h-full p-8">
      <div className="w-full max-w-xl space-y-6">
        {showVideo ? (
          <div className="space-y-4">
            <div className="rounded-2xl overflow-hidden border-2 border-border shadow-lg bg-card">
              <video
                key={videoKey}
                src="/mouse_open_field.mp4"
                // src="../../synthetic_mimic_mice_20251108_135756/mouse_open_field.mp4"
                controls 
                autoPlay
                className="w-full h-auto"
                poster="/video-placeholder.png"
              >
                Your browser does not support the video tag.
              </video>
            </div>
            <p className="text-sm text-muted-foreground text-center">
              Behavioral walk-through from the most recent simulation run.
            </p>
            <div className="flex justify-center">
              <Button variant="outline" onClick={() => setIsDialogOpen(true)}>
                View Prompting Instructions
              </Button>
            </div>
          </div>
        ) : (
          <div className="flex flex-col items-center space-y-4">
            <button
              type="button"
              onClick={() => setIsDialogOpen(true)}
              className="relative rounded-full focus:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2"
            >
              <div className="absolute inset-0 bg-primary/20 blur-2xl rounded-full" />
              <div className="relative h-48 w-48 bg-gradient-to-br from-primary/80 to-primary/40 rounded-full shadow-lg flex items-center justify-center">
                <div className="relative h-40 w-40">
                  <Image
                    src="/Brahma.png"
                    alt="Brahma AI logo"
                    fill
                    className="object-contain"
                    priority
                  />
                </div>
              </div>
            </button>
            <div className="text-center space-y-2">
              <h1 className="text-3xl font-bold text-foreground">Brahma AI</h1>
              <p className="text-muted-foreground text-sm">Tap the logo for prompting instructions</p>
              <p className="text-xs text-muted-foreground">
                Run a simulation to watch the behavioral walkthrough here.
              </p>
            </div>
          </div>
        )}
      </div>

      <Dialog open={isDialogOpen} onOpenChange={setIsDialogOpen}>
        <DialogContent className="max-w-2xl max-h-[80vh] overflow-y-auto">
          <DialogHeader>
            <DialogTitle className="text-2xl font-bold">How to Use Brahma AI</DialogTitle>
            <DialogDescription>
              Follow these steps to create and simulate synthetic mouse experiments
            </DialogDescription>
          </DialogHeader>
          
          <div className="space-y-4 mt-4">
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Dna className="h-5 w-5 text-primary" />
                  Step 1: Define Genotype
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-sm text-muted-foreground">
                  Specify the mouse strain and any genetic modifications. Example:
                </p>
                <pre className="mt-2 p-3 bg-muted rounded-md text-xs overflow-x-auto">
{`{
  "base": "C57BL/6J",
  "knockout": "Trp53"
}`}
                </pre>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Beaker className="h-5 w-5 text-primary" />
                  Step 2: Specify Treatment
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-sm text-muted-foreground">
                  Define experimental treatments or interventions. Example:
                </p>
                <pre className="mt-2 p-3 bg-muted rounded-md text-xs overflow-x-auto">
{`{
  "compound": "NewDrugX",
  "dose_mg_kg": 20,
  "days": 30
}`}
                </pre>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <BarChart3 className="h-5 w-5 text-primary" />
                  Step 3: Select Tests
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-sm text-muted-foreground">
                  Choose which tests to run on the synthetic mouse:
                </p>
                <ul className="mt-2 space-y-1 text-sm text-muted-foreground list-disc list-inside">
                  <li>Blood work analysis</li>
                  <li>Open field test (behavioral)</li>
                  <li>Elevated plus maze</li>
                  <li>Novel object recognition</li>
                </ul>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Mouse className="h-5 w-5 text-primary" />
                  Step 4: Run Simulation
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-sm text-muted-foreground">
                  Submit your request and receive a complete data package including:
                </p>
                <ul className="mt-2 space-y-1 text-sm text-muted-foreground list-disc list-inside">
                  <li>Genome sequences (FASTA)</li>
                  <li>Phenotype predictions (CSV)</li>
                  <li>Behavioral test results (CSV)</li>
                  <li>Biomarker data (JSON)</li>
                  <li>Summary reports</li>
                </ul>
              </CardContent>
            </Card>
          </div>

          <div className="mt-6 p-4 bg-muted/50 rounded-lg">
            <p className="text-xs text-muted-foreground">
              <strong>Tip:</strong> You can also specify a target mouse profile to mimic, 
              or let the system generate a realistic synthetic mouse based on your parameters.
            </p>
          </div>
        </DialogContent>
      </Dialog>
    </div>
  )
}

