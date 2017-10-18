macro import_fields(obj)
  quote
    eval(assign_fields($obj))
  end
end

function assign_fields(obj)
  expr = Expr(:block, [], Any)
  for field in fieldnames(obj)
    push!(expr.args, :($field = ($obj).$field))
  end
  return expr
end

export @import_fields
