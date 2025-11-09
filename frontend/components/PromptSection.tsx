"use client"

import { useMemo, useRef, useState } from "react"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import {
  Send,
  Loader2,
  Download,
  CheckCircle2,
  FileJson,
  Table as TableIcon,
  AlertCircle,
  Copy,
} from "lucide-react"

type SimulationResponse = {
  simulation_id: string
  download_url: string
  summary: {
    status: string
    phenotype_mimicry: boolean
    experiments_run: number
    tests_conducted: string[]
  }
  created_at: string
}

type UnpackedFile = {
  name: string
  type: "json" | "csv" | "text" | "other"
  content: string
  blob: Blob
  parsedJson?: unknown
  parsedCsv?: {
    headers: string[]
    rows: string[][]
  }
}

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000"

const renderJsonTree = (value: unknown, depth = 0): JSX.Element => {
  if (value === null || typeof value !== "object") {
    return (
      <span className="font-mono text-xs text-foreground/90">
        {typeof value === "string" ? value : JSON.stringify(value)}
      </span>
    )
  }

  if (Array.isArray(value)) {
    if (value.length === 0) {
      return (
        <div className="text-xs text-muted-foreground" style={{ marginLeft: depth * 12 }}>
          [empty array]
        </div>
      )
    }
    return (
      <div className="space-y-2" style={{ marginLeft: depth * 12 }}>
        {value.map((item, index) => (
          <div key={`${depth}-array-${index}`} className="border-l border-border/40 pl-3 space-y-1">
            <div className="text-[11px] font-medium text-muted-foreground">[{index}]</div>
            {renderJsonTree(item, depth + 1)}
          </div>
        ))}
      </div>
    )
  }

  const entries = Object.entries(value as Record<string, unknown>)
  if (!entries.length) {
    return (
      <div className="text-xs text-muted-foreground" style={{ marginLeft: depth * 12 }}>
        {'{ }'}
      </div>
    )
  }

  return (
    <div className="space-y-2" style={{ marginLeft: depth * 12 }}>
      {entries.map(([key, val]) => (
        <div key={`${depth}-obj-${key}`} className="border-l border-border/40 pl-3 space-y-1">
          <div className="text-xs font-semibold text-foreground">{key}</div>
          {renderJsonTree(val, depth + 1)}
        </div>
      ))}
    </div>
  )
}

export function PromptSection() {
  const [prompt, setPrompt] = useState("")
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [simulation, setSimulation] = useState<SimulationResponse | null>(null)
  const [unpackedFiles, setUnpackedFiles] = useState<UnpackedFile[]>([])
  const [zipBlob, setZipBlob] = useState<Blob | null>(null)
  const [isProcessingZip, setIsProcessingZip] = useState(false)
  const [jsonCopied, setJsonCopied] = useState<string | null>(null)
  const textareaRef = useRef<HTMLTextAreaElement>(null)

  const notifyVideoVisibility = (shouldShow: boolean) => {
    if (typeof window === "undefined") return
    window.dispatchEvent(new CustomEvent("simulation-video", { detail: { show: shouldShow } }))
  }

  const examplePrompts = useMemo(
    () => [
      {
        title: "Basic Mouse Simulation",
        prompt: `Create a C57BL/6J mouse with no modifications. Run blood work and open field test.`,
      },
      {
        title: "Drug Treatment Study",
        prompt: `Generate a mouse with Trp53 knockout. Apply NewDrugX at 20 mg/kg for 30 days. Run all tests.`,
      },
      {
        title: "Diet Intervention",
        prompt: `Create a mouse and simulate high-fat diet intervention. Analyze behavioral changes and biomarkers.`,
      },
    ],
    []
  )

  const determineFileType = (filename: string): UnpackedFile["type"] => {
    const lower = filename.toLowerCase()
    if (lower.endsWith(".json")) return "json"
    if (lower.endsWith(".csv")) return "csv"
    if (lower.endsWith(".txt") || lower.endsWith(".log") || lower.endsWith(".md")) return "text"
    return "other"
  }

  const parseCsv = (csvString: string) => {
    const lines = csvString.trim().split(/\r?\n/)
    if (!lines.length) return { headers: [] as string[], rows: [] as string[][] }
    const headers = lines[0].split(",").map((cell) => cell.trim())
    const rows = lines.slice(1).map((line) => line.split(",").map((cell) => cell.trim()))
    return { headers, rows }
  }

  const copyToClipboard = (text: string, id: string) => {
    if (typeof navigator === "undefined" || !navigator.clipboard) {
      console.warn("Clipboard API not available")
      return
    }
    navigator.clipboard.writeText(text).then(() => {
      setJsonCopied(id)
      setTimeout(() => setJsonCopied((prev) => (prev === id ? null : prev)), 2000)
    })
  }

  const fetchAndProcessZip = async (result: SimulationResponse) => {
    setIsProcessingZip(true)
    try {
      const downloadUrl = `${API_BASE_URL}${result.download_url}`
      const response = await fetch(downloadUrl)
      if (!response.ok) throw new Error("Failed to download simulation package.")
      const blob = await response.blob()
      setZipBlob(blob)

      const JSZip = (await import("jszip")).default
      const archive = await JSZip.loadAsync(blob)
      const files: UnpackedFile[] = []

      await Promise.all(
        Object.entries(archive.files).map(async ([filename, entry]) => {
          if (entry.dir) return
          const content = await entry.async("string")
          const type = determineFileType(filename)
          const file: UnpackedFile = {
            name: filename,
            type,
            content,
            blob: new Blob([content], {
              type:
                type === "json"
                  ? "application/json"
                  : type === "csv"
                  ? "text/csv"
                  : "text/plain",
            }),
          }

          if (type === "json") {
            try {
              file.parsedJson = JSON.parse(content)
            } catch {
              file.parsedJson = null
            }
          }

          if (type === "csv") {
            file.parsedCsv = parseCsv(content)
          }

          files.push(file)
        })
      )

      files.sort((a, b) => a.name.localeCompare(b.name))
      setUnpackedFiles(files)
    } finally {
      setIsProcessingZip(false)
    }
  }

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    if (!prompt.trim()) return
    setIsSubmitting(true)
    setError(null)
    setSimulation(null)
    setUnpackedFiles([])
    setZipBlob(null)
    notifyVideoVisibility(false)

    try {
      const lowerPrompt = prompt.toLowerCase()
      const tests: string[] = []
      if (lowerPrompt.includes("blood")) tests.push("blood_work")
      if (lowerPrompt.includes("open field") || lowerPrompt.includes("behavior")) tests.push("open_field")
      if (!tests.length) tests.push("blood_work", "open_field")

      let knockout: string | null = null
      if (lowerPrompt.includes("trp53")) knockout = "Trp53"

      const body = {
        genotype: { base: "C57BL/6J", knockout },
        treatment: {
          compound: lowerPrompt.includes("newdrugx") ? "NewDrugX" : "control",
          dose_mg_kg: lowerPrompt.includes("mg/kg") ? parseFloat(lowerPrompt.match(/(\d+(\.\d+)?)\s*mg\/kg/)?.[1] || "20") : 20,
          days: 30,
        },
        tests,
      }

      const res = await fetch(`${API_BASE_URL}/api/simulate`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      })

      if (!res.ok) {
        const data = await res.json().catch(() => ({}))
        throw new Error(data.detail || "Simulation failed. Please try again.")
      }

      const data: SimulationResponse = await res.json()
      setSimulation(data)
      await fetchAndProcessZip(data)
      notifyVideoVisibility(true)
    } catch (err) {
      console.error(err)
      setError(err instanceof Error ? err.message : "Unexpected error occurred.")
      notifyVideoVisibility(false)
    } finally {
      setIsSubmitting(false)
    }
  }

  const handleDownloadZip = async () => {
    if (zipBlob && simulation) {
      const url = URL.createObjectURL(zipBlob)
      const anchor = document.createElement("a")
      anchor.href = url
      anchor.download = `simulation_${simulation.simulation_id}.zip`
      document.body.appendChild(anchor)
      anchor.click()
      document.body.removeChild(anchor)
      URL.revokeObjectURL(url)
      return
    }
    if (simulation) {
      await fetchAndProcessZip(simulation)
    }
  }

  const handleExampleClick = (examplePrompt: string) => {
    setPrompt(examplePrompt)
    setError(null)
    setSimulation(null)
    setUnpackedFiles([])
    setZipBlob(null)
    notifyVideoVisibility(false)
    setTimeout(() => {
      textareaRef.current?.focus()
      textareaRef.current?.scrollIntoView({ behavior: "smooth", block: "center" })
    }, 80)
  }

  return (
    <div className="h-full flex flex-col p-8 space-y-6 overflow-auto">
      <Card className="flex-1 flex flex-col">
        <CardHeader>
          <CardTitle>Experiment Prompt</CardTitle>
          <CardDescription>
            Describe your experiment or tap an example. We’ll run the full virtual mouse simulation and unpack the results.
          </CardDescription>
        </CardHeader>
        <CardContent className="flex-1 flex flex-col space-y-4">
          <form onSubmit={handleSubmit} className="flex-1 flex flex-col space-y-4">
            <div className="flex-1">
              <textarea
                ref={textareaRef}
                value={prompt}
                onChange={(e) => setPrompt(e.target.value)}
                placeholder="Example: Create a C57BL/6J mouse with Trp53 knockout. Apply NewDrugX at 20 mg/kg for 30 days. Run blood work and open field test..."
                className="w-full h-full min-h-[220px] p-4 border-2 border-border rounded-md bg-background text-foreground resize-none focus:outline-none focus:ring-2 focus:ring-ring focus:ring-offset-2 focus:border-ring"
                disabled={isSubmitting}
              />
            </div>
            <div className="flex items-center justify-between">
              <Button type="submit" disabled={isSubmitting || !prompt.trim()} className="flex items-center gap-2">
                {isSubmitting ? (
                  <>
                    <Loader2 className="h-4 w-4 animate-spin" />
                    Running Simulation…
                  </>
                ) : (
                  <>
                    <Send className="h-4 w-4" />
                    Run Simulation
                  </>
                )}
              </Button>
            </div>
          </form>

          {error && (
            <Card className="border-destructive/40 bg-destructive/10">
              <CardContent className="p-4 flex gap-3 items-start">
                <AlertCircle className="h-5 w-5 text-destructive mt-0.5" />
                <div>
                  <p className="text-sm font-semibold text-destructive">Simulation failed</p>
                  <p className="text-sm text-destructive/80">{error}</p>
                </div>
              </CardContent>
            </Card>
          )}

          {simulation && (
            <div className="space-y-4">
              <Card className="border-border/70 bg-muted/40">
                <CardContent className="p-5 space-y-4">
                  <div className="flex items-start gap-3">
                    <CheckCircle2 className="h-5 w-5 text-primary mt-1" />
                    <div className="space-y-1 text-sm text-muted-foreground">
                      <p className="font-semibold text-foreground text-base">Simulation completed successfully!</p>
                      <p>Status: <span className="font-medium text-foreground">{simulation.summary.status}</span></p>
                      <p>Experiments run: <span className="font-medium text-foreground">{simulation.summary.experiments_run}</span></p>
                      <p>Tests conducted: <span className="font-medium text-foreground">{simulation.summary.tests_conducted.join(", ")}</span></p>
                      <p>Simulation ID: <span className="font-mono text-xs text-foreground/70">{simulation.simulation_id}</span></p>
                      <p className="text-xs">Completed at: {new Date(simulation.created_at).toLocaleString()}</p>
                    </div>
                  </div>
                  <div className="flex flex-col sm:flex-row gap-2">
                    <Button onClick={handleDownloadZip} disabled={isProcessingZip} className="flex-1 flex items-center justify-center gap-2">
                      {isProcessingZip ? (
                        <>
                          <Loader2 className="h-4 w-4 animate-spin" />
                          Preparing ZIP…
                        </>
                      ) : (
                        <>
                          <Download className="h-4 w-4" />
                          Download Results (ZIP)
                        </>
                      )}
                    </Button>
                  </div>
                </CardContent>
              </Card>

              {isProcessingZip && (
                <Card className="border-border/60 border-dashed bg-muted/30">
                  <CardContent className="p-4 flex items-center gap-3 text-sm text-muted-foreground">
                    <Loader2 className="h-4 w-4 animate-spin text-primary" />
                    Fetching and unpacking simulation package…
                  </CardContent>
                </Card>
              )}

              {!isProcessingZip && unpackedFiles.length > 0 && (
                <div className="space-y-3">
                  <h3 className="text-sm font-semibold text-foreground">Generated Data</h3>
                  {unpackedFiles.map((file) => {
                    const fileId = `file-${file.name}`
                    return (
                      <Card key={file.name} className="border-border/70">
                        <CardHeader className="flex flex-col gap-2 sm:flex-row sm:items-center sm:justify-between">
                          <div>
                            <CardTitle className="text-base font-semibold break-words">{file.name}</CardTitle>
                            <CardDescription className="capitalize text-xs">{file.type} file</CardDescription>
                          </div>
                          <div className="flex flex-wrap gap-2">
                            <Button variant="outline" size="sm" onClick={() => copyToClipboard(file.content, fileId)}>
                              <Copy className="mr-1 h-4 w-4" />
                              {jsonCopied === fileId ? "Copied!" : "Copy"}
                            </Button>
                            <Button
                              variant="outline"
                              size="sm"
                              onClick={() => {
                                const url = URL.createObjectURL(file.blob)
                                const anchor = document.createElement("a")
                                anchor.href = url
                                anchor.download = file.name.split("/").pop() || file.name
                                document.body.appendChild(anchor)
                                anchor.click()
                                document.body.removeChild(anchor)
                                URL.revokeObjectURL(url)
                              }}
                            >
                              <Download className="mr-1 h-4 w-4" />
                              Download
                            </Button>
                          </div>
                        </CardHeader>
                        <CardContent>
                          {file.type === "json" && (
                            <div className="space-y-2">
                              <div className="flex items-center gap-2 text-sm font-medium text-foreground">
                                <FileJson className="h-4 w-4" />
                                JSON Preview
                              </div>
                              {(() => {
                                const jsonData =
                                  file.parsedJson ??
                                  (() => {
                                    try {
                                      return JSON.parse(file.content)
                                    } catch {
                                      return null
                                    }
                                  })()
                                if (!jsonData) {
                                  return (
                                    <pre className="overflow-x-auto rounded-md border border-border bg-background/60 p-4 text-xs font-mono text-foreground whitespace-pre-wrap break-words">
                                      {file.content}
                                    </pre>
                                  )
                                }
                                return (
                                  <div className="rounded-md border border-border bg-background/60 p-4 text-sm space-y-2">
                                    {renderJsonTree(jsonData)}
                                  </div>
                                )
                              })()}
                            </div>
                          )}

                          {file.type === "csv" && file.parsedCsv && (
                            <div className="space-y-2">
                              <div className="flex items-center gap-2 text-sm font-medium text-foreground">
                                <TableIcon className="h-4 w-4" />
                                CSV Preview
                              </div>
                              <div className="overflow-x-auto rounded-md border border-border bg-background/60">
                                <table className="min-w-full divide-y divide-border text-xs">
                                  <thead className="bg-muted/60">
                                    <tr>
                                      {file.parsedCsv.headers.map((header, index) => (
                                        <th key={`${file.name}-header-${index}`} className="px-3 py-2 text-left font-semibold text-foreground">
                                          {header}
                                        </th>
                                      ))}
                                    </tr>
                                  </thead>
                                  <tbody className="divide-y divide-border">
                                    {file.parsedCsv.rows.length === 0 && (
                                      <tr>
                                        <td colSpan={file.parsedCsv.headers.length} className="px-3 py-2 text-muted-foreground">
                                          No data rows in this CSV file.
                                        </td>
                                      </tr>
                                    )}
                                    {file.parsedCsv.rows.map((row, rowIndex) => (
                                      <tr key={`${file.name}-row-${rowIndex}`} className="bg-background">
                                        {file.parsedCsv?.headers.map((_, colIndex) => (
                                          <td key={`${file.name}-cell-${rowIndex}-${colIndex}`} className="px-3 py-2">
                                            {row[colIndex] ?? ""}
                                          </td>
                                        ))}
                                      </tr>
                                    ))}
                                  </tbody>
                                </table>
                              </div>
                            </div>
                          )}

                          {(file.type === "text" || file.type === "other") && (
                            <pre className="overflow-x-auto rounded-md border border-border bg-background/60 p-4 text-xs font-mono text-foreground whitespace-pre-wrap break-words">
                              {file.content.length > 5000 ? `${file.content.slice(0, 5000)}…` : file.content}
                            </pre>
                          )}
                        </CardContent>
                      </Card>
                    )
                  })}
                </div>
              )}
            </div>
          )}

          <div className="space-y-3">
            <h3 className="text-sm font-semibold text-foreground">Example Prompts</h3>
            <div className="space-y-2">
              {examplePrompts.map((example, index) => (
                <Card
                  key={index}
                  className="cursor-pointer hover:bg-accent/50 transition-colors active:bg-accent/70"
                  onClick={() => handleExampleClick(example.prompt)}
                >
                  <CardContent className="p-4">
                    <p className="text-sm font-medium mb-1">{example.title}</p>
                    <p className="text-xs text-muted-foreground">{example.prompt}</p>
                  </CardContent>
                </Card>
              ))}
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}

