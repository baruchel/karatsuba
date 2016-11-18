def compile_plan(name, d):
    context = {}
    exec (d%name) in context
    return context[name]
