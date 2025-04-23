using TestItems

@testitem "aqua" begin
  using UnderwaterAcoustics
  using Aqua
  Aqua.test_all(UnderwaterAcoustics; persistent_tasks=false)
end
