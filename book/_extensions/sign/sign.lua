function sign(args)
  local var = pandoc.utils.stringify(args[1])
  return pandoc.Str("See @" .. var)
end