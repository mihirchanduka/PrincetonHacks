import { LogoSection } from "@/components/LogoSection"
import { PromptSection } from "@/components/PromptSection"

export default function Home() {
  return (
    <div className="min-h-screen flex">
      {/* Left Side - Logo and Instructions */}
      <div className="w-1/2 border-r border-border bg-gradient-to-br from-background to-muted/20">
        <LogoSection />
      </div>

      {/* Right Side - Prompt Area */}
      <div className="w-1/2 bg-background">
        <PromptSection />
      </div>
    </div>
  )
}


